%==========================================================================
% This script analyzes the temporal stability of egocentric bearing cells
% (ELRPD cells).
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.myRNG             = 333;
param.clusterName       = 'HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768; % maximum absolute yaw value in the arena
param.newTimeRes        = 0.1; % temporal resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates   	= 101; % number of surrogates
%- for data parts
param.numDataParts      = 2; % median split into early and late parts
%- for direction
direc                   = [];
direc.angularRes        = 30; % angular resolution
direc.angleEdges        = deg2rad(-180:direc.angularRes:180);
direc.angleCenters      = transpose(movmean(direc.angleEdges, 2, 'endpoints', 'discard'));
%- for location
place                   = [];
place.locRes            = 10; % spatial resolution
place.idxTemplate       = flipud(reshape(1:place.locRes^2, place.locRes, place.locRes));
place.xEdges            = linspace(-4500, 4500, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.yEdges            = linspace(-4500, 4500, place.locRes + 1);
place.yCenters          = movmean(place.yEdges, 2, 'endpoints', 'discard');
place.cutoffLocBin      = 0.95;
%- for egocentric local directions
locDir                  = [];
locDir.locRes           = 12; % spatial resolution
locDir.xEdges           = linspace(-5400, 5400, locDir.locRes + 1);
locDir.xCenters         = movmean(locDir.xEdges, 2, 'endpoints', 'discard');
locDir.yEdges           = linspace(-5400, 5400, locDir.locRes + 1);
locDir.yCenters         = movmean(locDir.yEdges, 2, 'endpoints', 'discard');
locDir.angularRes       = 30; % angular resolution
locDir.angleEdges       = deg2rad(-180:locDir.angularRes:180);
locDir.angleCenters     = transpose(movmean(locDir.angleEdges, 2, 'endpoints', 'discard'));
%- local reference points
tmpX                    = meshgrid(locDir.xCenters);
tmpY                    = transpose(meshgrid(locDir.yCenters));
locDir.LRPs             = [tmpX(:), tmpY(:)];
% exclude LRPs that are too far from center
D_LRPs2Ctr              = sqrt(locDir.LRPs(:, 1) .^ 2 + locDir.LRPs(:, 2) .^ 2);
locDir.LRPs2Exclude     = D_LRPs2Ctr > 5200;
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sum of squares

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.save      = strcat('E:\OpenField\ELRPDCellTemporalStability_20200727\20200904_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
    '_EgoLocDir', num2str(locDir.angularRes), '\', ...
    'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotation_Smooth', num2str(param.naviSmooth), '\');
mkdir(paths.save);

% subjects
subjects    = {...
    'OF_001a'; ...
    'OF_001b'; ...
    'OF_002'; ...
    'OF_003a'; ...
    'OF_003b'; ...
    'OF_004a'; ...
    'OF_004b'; ...
    'OF_005'; ...
    'OF_006'; ...
    'OF_007'; ...
    'OF_008'; ...
    'OF_009'; ...
    'OF_010'; ...
    'OF_011'; ...
    'OF_012'; ...
    'OF_013a'; ...
    'OF_013b'; ...
    'OF_014'; ...
    };

%% save settings
save(strcat(paths.save, 'settings'));

%% load previous ELRPD cell analysis

% results from ELRPD cell analysis
ELRPD  = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes');

%% preallocations

% main results for each cell
allRes    	= [];

% bookkeeping
exByWC    	= []; % excluded by wave-clus
exByRefAna	= []; % excluded by reference analysis

% for behavioral control analyses
allBehLRP 	= [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(param.myRNG);
        
    %% wires to investigate
    
    % available microwires
    wires   = dir(fullfile(paths.spike, subjects{iSub}, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    %% loop through wires
    for iWire = 1:size(wires, 1)
        
        % report
        fprintf('\tWire: %s.\n', wires(iWire).name);
        
        %% brain region and MNI coordinates
        
        % subject information
        s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion  = s.subjectdata.micro2macro{logIdx, 3}; % brain region
        wireMNI     = s.subjectdata.micro2macro{logIdx, 5}; % MNI coordinates
        
        %% spike times in behavioral time
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
        catch
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC  = cat(1, exByWC, [iSub, iWire]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on inspection)
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', param.clusterName, '.mat'));
        
        % load behavioral times of spikes
        ccb = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('Size of "cluster_class" not congruent with "cluster_class_behtime".');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
                        
            % continue if this cluster was not included in the ELRPD-cell
            % analysis or if it is not an ELRPD cell
            refIdx  = all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2);
            if sum(refIdx) ~= 1 || ~(ELRPD.allRes(refIdx).locdir_clusRank > 0.95)
                fprintf('\t\t- This cluster was not included in the ELRPD-cell analysis or is not an ELRPD cell.\n');
                exByRefAna  = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% process behavioral data
                        
            % this cell's reference point
            thisCOMxy           = ELRPD.allRes(refIdx).locdir_COMxy;
            
            % process behavior with respect to local reference point
            dt                  = [];
            dt.paths            = paths;
            dt.subject          = subjects{iSub};
            dt.direc            = direc;
            dt.place            = place;
            dt.locDir           = locDir;
            dt.locDir.refPoint  = thisCOMxy; % specific LRP
            dt.maxAbsYaw        = param.maxAbsYaw;
            dt.maxNumTrials     = param.maxNumTrials;
            dt.newTimeRes       = param.newTimeRes;
            dt.naviCutoff       = param.naviCutoff;
            dt.naviSmooth       = param.naviSmooth;
            dt.rotSmooth        = param.rotSmooth;
            tmpOut              = LK_EvalBehOverall_20200623(dt);
            
            % use new timeline for the analyses
            outBehLRP           = tmpOut.new;
            
            % collect behavioral output for all subjects
            allBehLRP           = cat(1, allBehLRP, outBehLRP');
            
            %% mask for analysis timepoints (timepoints to be included in the analysis)
            
            % boolean that encodes whether specific time bins shall be
            % included in the analysis
            bMask4Analysis      = outBehLRP.bNavi | outBehLRP.bRotation;
            fprintf('Information: %.3f%% of the "active time" is included in the analysis.\n', ...
                sum(outBehLRP.durations(bMask4Analysis)) / sum(outBehLRP.durations(outBehLRP.bActive)) * 100);
            
            %% number of spikes and FR per time bin
                        
            % number of spikes per timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBehLRP.time), 0]); % add 0 spikes for the last time-bin
            
            % FR per timepoint
            FR          = numSpikes ./ outBehLRP.durations;
            
            %% ANOVA for modulation of firing rate by egocentric direction towards local reference points, entire data
            
            % compute ANOVAN
            dt              = [];
            dt.FR           = FR(bMask4Analysis);
            dt.group        = [outBehLRP.yawBin(bMask4Analysis), outBehLRP.xyBin(bMask4Analysis), outBehLRP.egoLocDirBin(bMask4Analysis)];
            dt.groupNames   = {'D'; 'P'; 'ELRPD'}; % predictors
            dt.ANOVA        = ANOVA;
            empANOVA        = LK_CompANOVA_070719(dt);
                        
            % add additional mask to only analyze group levels that are
            % also included in the main ANOVA
            bAddMask4Analysis       = ismember(outBehLRP.yawBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')})) & ...
                ismember(outBehLRP.xyBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'P')})) & ...
                ismember(outBehLRP.egoLocDirBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'ELRPD')}));            
            
            % update analysis mask
            bMask4Analysis  = bMask4Analysis & bAddMask4Analysis;
            
            %% index for different data parts
            
            % divide data into parts
            cs      = cumsum(bMask4Analysis) ./ max(cumsum(bMask4Analysis));
            idxPart = ones(size(cs, 1), 1);
            for iPart = 1:param.numDataParts - 1
                idxPart(cs > iPart / param.numDataParts)    = iPart + 1;
            end
                                                
            %% ANOVA for modulation of firing rate by egocentric direction towards local reference points, separately for the different data parts
            
            % for each data part, estimate modulation of firing rate by
            % egocentric direction towards the reference point
            locdir_F        = nan(param.numDataParts, 1); % size = 2 x 1
            locdir_corrFR   = nan(param.numDataParts, numel(locDir.angleCenters)); % size = 2 x 12
            for iPart = 1:param.numDataParts
                
                % select data points from this data part
                thisBMask4Analysis          = bMask4Analysis & idxPart == iPart;
                
                % compute ANOVAN
                dt                          = [];
                dt.FR                       = FR(thisBMask4Analysis);
                dt.group                    = [outBehLRP.yawBin(thisBMask4Analysis), outBehLRP.xyBin(thisBMask4Analysis), outBehLRP.egoLocDirBin(thisBMask4Analysis)];
                dt.groupNames               = {'D'; 'P'; 'ELRPD'};
                dt.ANOVA                    = ANOVA;
                dt.ANOVA.minNumObsPerBin    = 1; % no additional requirement as compared to the overall ANOVA
                dt.bSkipA1                  = 'yes';
                empANOVA                    = LK_CompANOVA_070719(dt);
                
                % get output F value
                locdir_F(iPart, 1)          = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'ELRPD'), strcmp(empANOVA.tbl(1, :), 'F')));
                
                % corrected firing rate for ELRPD
                tmpFR                       = empANOVA.m{strcmp(empANOVA.groupNames, 'ELRPD')}(:, 1);
                idx                         = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'ELRPD')});
                locdir_corrFR(iPart, idx)   = tmpFR;
            end
            
            %% collect information for this unit
            
            % basics
            unitRes                 = [];
            unitRes.idx             = [iSub, iWire, iClus];
            unitRes.bMask4Analysis  = bMask4Analysis;
            unitRes.idxPart         = idxPart;
            unitRes.wireRegion      = wireRegion;
            unitRes.wireMNI         = wireMNI;
            
            % ELRPD
            unitRes.COMxy           = thisCOMxy;
            unitRes.locdir_F        = locdir_F;
            unitRes.locdir_corrFR   = locdir_corrFR;
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);
            
            %% close all open figures
            close all;
            toc            
        end
    end
end

%% save all results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));

% all cells' indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');
fprintf('Number of cells: %d.\n', size(allUnitIdx, 1));

%% temporal stability of tuning curves

% loop through cells
prefELRPD   = nan(size(r.allRes, 1), r.param.numDataParts); % preferred ELRPD in both data halves
tempStab    = nan(size(r.allRes, 1), 1); % temporal stability
for iCell = 1:size(r.allRes, 1)

    % correlation between tuning curves
    rho                 = corr(transpose(r.allRes(iCell).locdir_corrFR));
    logMat              = tril(ones(size(rho)), -1) ~= 0; % lower triangle
    tempStab(iCell, 1)  = mean(rho(logMat));
    
    % preferred ELRPD in all data parts
    for iPart = 1:r.param.numDataParts
        prefELRPD(iCell, iPart) = circ_mean(r.locDir.angleCenters, transpose(r.allRes(iCell).locdir_corrFR(iPart, :)));
    end
end

% circular-circular correlation between preferred ELRPDs
[rho, pval] = circ_corrcc(prefELRPD(:, 1), prefELRPD(:, 2));
fprintf('Circular-circular correlation between preferred egocentric bearings in both data halves: rho = %.3f, p = %.3f.\n', rho, pval);

% t-test of correlation values against zero
[~, p, ~, stats]    = ttest(tempStab);
fprintf('Temporal stability of ELRPD cells: mean r = %.3f, SEM = %.3f, t(%d) = %.3f, p = %.3f.\n', ...
    mean(tempStab), std(tempStab) / sqrt(numel(tempStab)), stats.df, stats.tstat, p);

% figure
f = figure('units', 'centimeters', 'position', [2, 22, 8, 8]);
hold on;
histogram(tempStab, -1:0.25:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
tmpAx = get(gca);
plot([mean(tempStab), mean(tempStab)], [0, max(tmpAx.YLim)], ...
    'Color', [0, 0, 0], 'LineWidth', 2);
hold off;
xl = xlabel('Temporal stability (\itr\rm)');
yl = ylabel('Count');
set(gca, ...
    'xtick', -1:0.5:1, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ELRPDCells_temporalStability_20200723'), '-dtiff', '-r450');
