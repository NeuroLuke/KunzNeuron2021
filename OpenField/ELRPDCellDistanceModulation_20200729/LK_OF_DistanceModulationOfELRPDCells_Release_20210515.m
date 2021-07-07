%==========================================================================
% This script analyzes whether egocentric bearing cells (ELRPD cells)
% exhibit distance modulation towards their reference point (linear
% distance tuning and conjunctive bearing-distance tuning).
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% variables
param                   = [];
param.myRNG             = 9999;
param.clusterName       = 'HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768; % maximum absolute yaw value
param.newTimeRes        = 0.1; % new time resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates   	= 101; % number of surrogates
param.tailAlpha         = 0.025; % alpha-level at both ends of the surrogate distribution
param.tailThresh        = 1 - param.tailAlpha; % rank-threshold at both ends of the surrogate distribution
%- for direction
direc                   = [];
direc.angularRes        = 30; % angular resolution
direc.angleEdges        = deg2rad(-180:direc.angularRes:180);
direc.angleCenters      = transpose(movmean(direc.angleEdges, 2, 'endpoints', 'discard'));
%- for location
place                   = [];
place.locRes            = 10; % spatial resolution
place.idxTemplate       = flipud(reshape(1:place.locRes^2, place.locRes, place.locRes)); % template for location indexing
place.xEdges            = linspace(-4500, 4500, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.yEdges            = linspace(-4500, 4500, place.locRes + 1);
place.yCenters          = movmean(place.yEdges, 2, 'endpoints', 'discard');
place.cutoffLocBin      = 0.95;
%- for local directions
locDir                  = [];
locDir.locRes           = 12; % spatial resolution
locDir.xEdges           = linspace(-5400, 5400, locDir.locRes + 1);
locDir.xCenters         = movmean(locDir.xEdges, 2, 'endpoints', 'discard');
locDir.yEdges           = linspace(-5400, 5400, locDir.locRes + 1);
locDir.yCenters         = movmean(locDir.yEdges, 2, 'endpoints', 'discard');
locDir.angularRes       = 30; % angular resolution
locDir.angleEdges       = deg2rad(-180:locDir.angularRes:180);
locDir.angleCenters     = transpose(movmean(locDir.angleEdges, 2, 'endpoints', 'discard'));
%- for distance
distnc                  = [];
distnc.corrType         = 'Pearson'; % correlation type
distnc.res              = 500; % extent of bins (for graphical depiction)
distnc.binEdges         = 0:distnc.res:10000;
distnc.binCenters       = movmean(distnc.binEdges, 2, 'endpoints', 'discard');
%- for bearing-distance fields
bearDist                    = [];
bearDist.numSurrogates  	= 1001; % number of surrogates
bearDist.distanceBinEdges 	= 0:200:10000; % 50 distance bins
bearDist.distanceBinCenters = movmean(bearDist.distanceBinEdges, 2, 'endpoints', 'discard');
bearDist.bearingBinEdges    = linspace(-pi, pi, 25); % 24 bearing bins
bearDist.bearingBinCenters  = movmean(bearDist.bearingBinEdges, 2, 'endpoints', 'discard');
bearDist.smoothFac          = 5;
bearDist.smoothType         = 'gaussian';
%- for ANOVA procedure
ANOVA                       = [];
ANOVA.minNumObsPerBin       = 5; % minimum number of observations per bin
ANOVA.model                 = 'linear'; % type of ANOVA model
ANOVA.sstype                = 2; % type of sums of squares

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.arena     = 'E:\OpenField\Arena\';
paths.save      = strcat('E:\OpenField\ELRPDCellDistanceModulation_20200729\20201013_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
    '_EgoLocDir', num2str(locDir.angularRes), '_Dist', num2str(distnc.res), '\', ...
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

%% reference analysis

% load reference analysis (to only include ELRPD cells in the analysis)
refAna  = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'locDir');

%% preallocations

% main results for each cell
allRes      = [];

% for sanity checking
exByWC      = []; % excluded by wave-clus
exByRefAna  = []; % excluded by reference analysis

% all behavioral data
allBehLRP   = [];

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
        
        %% brain region of this wire
        
        % subject information
        s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion  = s.subjectdata.micro2macro{logIdx, 3};
        
        %% get spike times in behavioral time
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
        catch
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC  = cat(1, exByWC, [iSub, iWire]); % bookkeeping
            continue;
        end
        
        % load behavioral times of spikes
        ccb     = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('Size of "cluster_class" not congruent with "cluster_class_behtime".');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
            
            % continue if this cluster was not included in the reference
            % analysis or if it is not an ELRPD cell
            refIdx  = all([iSub, iWire, iClus] == cell2mat({refAna.allRes.idx}'), 2);
            if sum(refIdx) ~= 1 || ~(refAna.allRes(refIdx).locdir_clusRank > 0.95)
                fprintf('\t\t- This cluster was not included in the reference analysis or is not an ELPRD cell.\n');
                exByRefAna  = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % cluster-specific save path
            clusSavePath    = strcat(paths.save, subjects{iSub}, '_chan', num2str(iWire), '_Clus', num2str(iClus), '\');
            if ~exist(clusSavePath, 'dir')
                mkdir(clusSavePath);
            end
            
            % data for this cluster
            thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% behavioral data for this cell's LRP
            
            % this cell's reference point
            locdir_COMxy    = refAna.allRes(refIdx).locdir_COMxy;
            
            % this cell's preferred ELRPD bearing
            LRPIdx          = all(refAna.allRes(refIdx).locdir_COMxyClosestLRP == refAna.locDir.LRPs, 2);
            prefELRPD       = circ_mean(refAna.locDir.angleCenters, transpose(refAna.allRes(refIdx).locdir_corrFR(LRPIdx, :)));
            
            % process behavior with respect to local reference point
            behLRPFile      = strcat(clusSavePath, 'outBehLRP.mat');
            if exist(behLRPFile, 'file')
                
                % load previously saved results
                fprintf('Loading behavioral data from "%s".\n', behLRPFile);
                tmp      	= load(behLRPFile);
                outBehLRP  	= tmp.outBehLRP;
            else
                
                % process behavioral information
                dt                  = [];
                dt.paths            = paths;
                dt.subject          = subjects{iSub};
                dt.direc            = direc;
                dt.place            = place;
                dt.locDir           = locDir;
                dt.locDir.refPoint  = locdir_COMxy; % reference point
                dt.maxAbsYaw        = param.maxAbsYaw;
                dt.maxNumTrials     = param.maxNumTrials;
                dt.newTimeRes       = param.newTimeRes;
                dt.naviCutoff       = param.naviCutoff;
                dt.naviSmooth       = param.naviSmooth;
                dt.rotSmooth        = param.rotSmooth;
                tmpOut              = LK_EvalBehOverall_20200623(dt);
                
                % use behavioral data with new timeline
                outBehLRP           = tmpOut.new;
                
                % save output
                save(behLRPFile, 'outBehLRP');
            end
            
            % collect all behavioral information
            allBehLRP   = cat(1, allBehLRP, outBehLRP');
            
            %% calculate Euclidean distance between subject's location and COM
            
            % for each timepoint, distance between current location and the
            % reference point
            dist2COM   	= pdist2(outBehLRP.behinfo(:, 2:3), locdir_COMxy);
            
            % discretize the distances for use in the ANOVA
            dist2COMBin = discretize(dist2COM, distnc.binEdges);
            
            %% mask for analysis timepoints
            
            % boolean that encodes whether specific time bins shall be
            % included in the analysis
            bMask4Analysis  = outBehLRP.bNavi | outBehLRP.bRotation;
            fprintf('Information: %.3f%% of the "active time" is included into the analysis.\n', ...
                sum(outBehLRP.durations(bMask4Analysis)) / sum(outBehLRP.durations(outBehLRP.bActive)) * 100);
            
            % sanity check
            if any(bMask4Analysis ~= refAna.allRes(refIdx).bMask4Analysis)
                error('Current "bMask4Analysis" is not congruent with "bMask4Analysis" from the reference analysis.\n');
            end
            
            %% number of spikes and FR per time bin
            
            % number of spikes and FR for each behavioral timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBehLRP.time), 0]); % add 0 spikes for the last time-bin
            FR          = numSpikes ./ outBehLRP.durations;
            
            %% empirical ANOVA to estimate FR as a function of distance to the reference point
            
            % for each local reference point, estimate modulation of firing
            % rate by distance towards this LRP
            resFile = strcat(clusSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_locdir_emp.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                     = load(resFile);
                dist2COM_F              = tmp.dist2COM_F; % F-value
                dist2COM_corrFR         = tmp.dist2COM_corrFR; % tuning curve
                
                % further data
                residEmp                = tmp.residEmp; % empirical residuals
                bMask4Analysis          = tmp.bMask4Analysis;
                bAddMask4Analysis       = tmp.bAddMask4Analysis;
            else
                
                %--- compute 4-way ANOVA with distance as a factor (for
                % graphical depiction)
                dt              = [];
                dt.FR           = FR(bMask4Analysis);
                dt.group        = [outBehLRP.yawBin(bMask4Analysis), outBehLRP.xyBin(bMask4Analysis), outBehLRP.egoLocDirBin(bMask4Analysis), ...
                    dist2COMBin(bMask4Analysis)];
                dt.groupNames   = {'D'; 'P'; 'ELRPD'; 'dist2COM'}; % predictors
                dt.ANOVA        = ANOVA;
                dt.bSkipA1      = 'yes';
                empANOVA        = LK_CompANOVA_070719(dt);
                
                % output F value
                dist2COM_F            	= cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'dist2COM'), strcmp(empANOVA.tbl(1, :), 'F')));
                
                % tuning curve for dist2COM (for plotting)
                dist2COM_corrFR       	= nan(1, numel(distnc.binCenters));
                tmpFR                 	= empANOVA.m{strcmp(empANOVA.groupNames, 'dist2COM')}(:, 1);
                idx                  	= cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'dist2COM')});
                dist2COM_corrFR(idx)  	= tmpFR;
                
                %--- obtain residuals from 3-way ANOVA for later use
                dt              = [];
                dt.FR           = FR(bMask4Analysis);
                dt.group        = [outBehLRP.yawBin(bMask4Analysis), outBehLRP.xyBin(bMask4Analysis), outBehLRP.egoLocDirBin(bMask4Analysis)];
                dt.groupNames   = {'D'; 'P'; 'ELRPD'};
                dt.ANOVA        = ANOVA;
                dt.bSkipA1      = 'yes';
                empANOVA        = LK_CompANOVA_070719(dt);
                residEmp     	= empANOVA.stats.resid; % residuals after controlling for direction, place, and ELRPD
                
                % additional mask to analyze the same samples as in the
                % ANOVA
                bAddMask4Analysis   = ismember(outBehLRP.yawBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')})) & ...
                    ismember(outBehLRP.xyBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'P')})) & ...
                    ismember(outBehLRP.egoLocDirBin, cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'ELRPD')}));
                
                % save result
                save(resFile, 'dist2COM_F', 'dist2COM_corrFR', 'residEmp', 'bMask4Analysis', 'bAddMask4Analysis');
            end
            
            %% surrogates
            
            % random numbers for circular shift of firing rates
            randShifts  = transpose(datasample(1:numel(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
            
            % create surrogates or load previously saved surrogates
            resFile     = strcat(clusSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                 	= load(resFile);
                dist2COM_Fsurro      	= tmp.dist2COM_Fsurro;
                dist2COM_corrFRsurro    = tmp.dist2COM_corrFRsurro;
                residSurro              = tmp.residSurro;
            else
                
                % preallocate
                dist2COM_Fsurro       	= nan(param.numSurrogates, 1);
                dist2COM_corrFRsurro  	= nan(param.numSurrogates, numel(distnc.binCenters));
                residSurro              = nan(size(residEmp, 1), param.numSurrogates);
                
                % loop through surrogates
                parfor iSurro = 1:param.numSurrogates
                    
                    % create surrogate firing rate by circularly shifting
                    % firing rates relative to the behavioral data
                    FRsurro         = circshift(FR(bMask4Analysis), randShifts(iSurro, 1)); % mask firing rates by analysis time window
                    
                    %% surrogate ANOVA for modulation of firing rate by egocentric direction towards local reference points
                    
                    %--- compute 4-way ANOVAN with distance as a factor
                    dt              = [];
                    dt.FR           = FRsurro; % you already masked the FR by the analysis period (see above)
                    dt.group        = [outBehLRP.yawBin(bMask4Analysis), outBehLRP.xyBin(bMask4Analysis), outBehLRP.egoLocDirBin(bMask4Analysis), ...
                        dist2COMBin(bMask4Analysis)];
                    dt.groupNames   = {'D'; 'P'; 'ELRPD'; 'dist2COM'};
                    dt.ANOVA        = ANOVA;
                    dt.bSkipA1      = 'yes';
                    surroANOVA      = LK_CompANOVA_070719(dt);
                    
                    % output F value
                    dist2COM_Fsurro(iSurro, 1)  = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'dist2COM'), strcmp(surroANOVA.tbl(1, :), 'F')));
                    
                    % tuning curve for dist2COM
                    tmp_dist2COM_corrFRsurro    	= nan(1, numel(distnc.binCenters));
                    tmpFR                           = surroANOVA.m{strcmp(surroANOVA.groupNames, 'dist2COM')}(:, 1);
                    idx                             = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'dist2COM')});
                    tmp_dist2COM_corrFRsurro(idx)   = tmpFR;
                    
                    % collect across surrogates
                    dist2COM_corrFRsurro(iSurro, :) = tmp_dist2COM_corrFRsurro;
                    
                    %--- obtain residuals from 3-way ANOVA for later use
                    dt                      = [];
                    dt.FR                   = FRsurro; % you already masked the FR by the analysis period (see above)
                    dt.group                = [outBehLRP.yawBin(bMask4Analysis), outBehLRP.xyBin(bMask4Analysis), outBehLRP.egoLocDirBin(bMask4Analysis)];
                    dt.groupNames           = {'D'; 'P'; 'ELRPD'};
                    dt.ANOVA                = ANOVA;
                    dt.bSkipA1              = 'yes';
                    surroANOVA              = LK_CompANOVA_070719(dt);
                    residSurro(:, iSurro)   = surroANOVA.stats.resid; % residuals after controlling for direction, place, and ELRPD                                        
                end
                
                % save surrogates
                save(resFile, ...
                    'dist2COM_Fsurro', 'dist2COM_corrFRsurro', 'residSurro');
            end
            
            %% collect information for this unit
            
            % basics
            unitRes                         = [];
            unitRes.idx                     = [iSub, iWire, iClus];
            unitRes.thisSpike               = thisSpike;
            unitRes.t                       = t;
            unitRes.numSpikes               = numSpikes;
            unitRes.wireRegion              = wireRegion;
            
            % reference point and preferred ELRPD bearing
            unitRes.locdir_COMxy            = locdir_COMxy;
            unitRes.prefELRPD               = prefELRPD;
            
            % assessment of distance tuning from ANOVA
            unitRes.dist2COM_F              = dist2COM_F;
            unitRes.dist2COM_Fsurro         = dist2COM_Fsurro;
            unitRes.dist2COM_corrFR         = dist2COM_corrFR;
            unitRes.dist2COM_corrFRsurro    = dist2COM_corrFRsurro;
            
            % relevant data
            unitRes.bMask4Analysis          = bMask4Analysis;
            unitRes.FR                      = FR;
            unitRes.dist2COM                = dist2COM;
            unitRes.bAddMask4Analysis       = bAddMask4Analysis;
            unitRes.residEmp                = residEmp; % empirical residuals
            unitRes.residSurro              = residSurro; % surrogate residuals
            
            % save unit result
            save(strcat(clusSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_unitRes'), 'unitRes');
            
            % collect across units
            allRes  = cat(1, allRes, unitRes);
            
            %% timing
            toc
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));

% all cell indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');
fprintf('Number of cells: %d.\n', size(allUnitIdx, 1));

%% identify distance-modulated ELRPD cells based on the correlations between all distance values and all residuals

% report
fprintf('\nEvaluation of distance-modulated ELRPD cells based on the correlation between distances and residuals from the 3way ANOVA.\n');

% preallocate and loop through cells
allRhoResid2Dist        = nan(size(r.allRes, 1), 1);
allRhoResid2DistSurro   = nan(size(r.allRes, 1), r.param.numSurrogates);
for iCell = 1:size(r.allRes, 1)
    
    % empirical correlation between distances and residuals from 3way ANOVA
    bTmpMask                    = r.allRes(iCell).bMask4Analysis & r.allRes(iCell).bAddMask4Analysis;
    allRhoResid2Dist(iCell, 1)  = corr(r.allRes(iCell).dist2COM(bTmpMask), r.allRes(iCell).residEmp, 'type', r.distnc.corrType);
    
    % surrogate correlation
    allRhoResid2DistSurro(iCell, :) = corr(r.allRes(iCell).dist2COM(bTmpMask), r.allRes(iCell).residSurro, 'type', r.distnc.corrType);
end

% positive distance cells
posRanksResid     	= sum(allRhoResid2Dist > allRhoResid2DistSurro, 2) ./ sum(~isnan(allRhoResid2DistSurro), 2);
bPosDistCellResid  	= posRanksResid > r.param.tailThresh;
fprintf('Number of positive distance cells: %d (P = %.3f).\n', sum(bPosDistCellResid), myBinomTest(sum(bPosDistCellResid), numel(bPosDistCellResid), r.param.tailAlpha));

% negative distance cells
negRanksResid     	= sum(allRhoResid2Dist < allRhoResid2DistSurro, 2) ./ sum(~isnan(allRhoResid2DistSurro), 2);
bNegDistCellResid 	= negRanksResid > r.param.tailThresh;
fprintf('Number of negative distance cells: %d (P = %.3f).\n', sum(bNegDistCellResid), myBinomTest(sum(bNegDistCellResid), numel(bNegDistCellResid), r.param.tailAlpha));

% combine positive and negative distance cells
bDistCellResid    	= bPosDistCellResid | bNegDistCellResid;
fprintf('Number of positive/negative distance cells: %d (P = %.3f).\n', sum(bDistCellResid), myBinomTest(sum(bDistCellResid), numel(bDistCellResid), r.param.tailAlpha * 2));

%% create figures for dependency of FRs on reference-point distance (using binned data for graphical depiction)

% preallocate linear fit
allLinFit   = cell(size(r.allRes, 1), 1);

% loop through cells
for iCell = 1:size(r.allRes, 1)
    
    %% create figure
    
    % convert rank into p-value
    if posRanksResid(iCell) > negRanksResid(iCell) % positive distance effect?
        pval4plot   = 1 - posRanksResid(iCell);
    elseif posRanksResid(iCell) < negRanksResid(iCell) % negative distance effect?
        pval4plot   = 1 - negRanksResid(iCell);
    end
    
    % create figure
    ld                  = [];
    ld.visible          = 'off';
    ld.FR               = r.allRes(iCell).dist2COM_corrFR;
    ld.distnc           = r.distnc;
    ld.locdir_COMxy     = r.allRes(iCell).locdir_COMxy;
    ld.prefELRPD        = r.allRes(iCell).prefELRPD;
    ld.thisSpike        = r.allRes(iCell).thisSpike;
    ld.sr               = r.allRes(iCell).t.par.sr;
    ld.nspk             = sum(r.allRes(iCell).numSpikes(r.allRes(iCell).bMask4Analysis));
    if pval4plot < 0.01
        ld.figTitle     = {'Distance tuning', ['\itP\rm < 0.01 (', r.allRes(iCell).wireRegion, ')']};
    elseif pval4plot < param.tailAlpha && round(pval4plot, 2) == param.tailAlpha
        ld.figTitle     = {'Distance tuning', ['\itP\rm < ', num2str(param.tailAlpha), ' (', r.allRes(iCell).wireRegion, ')']};
    else
        ld.figTitle     = {'Distance tuning', ['\itP\rm = ', num2str(pval4plot, '%.2f'), ' (', r.allRes(iCell).wireRegion, ')']};
    end
    [f, allLinFit{iCell}]   = LK_PlotLinearDistModulation_20201011(ld);
    
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(paths.save, r.subjects{r.allRes(iCell).idx(1)}, '_chan', num2str(r.allRes(iCell).idx(2)), ...
        '_Clus', num2str(r.allRes(iCell).idx(3)), '_Dist2COM'), '-dtiff', '-r600');
    close(f);
end

% linear fits of all units
allLinFit   = cell2mat(allLinFit);

%% summary plot on linear distance tuning

% normalize firing rates (so that max = 1 and min = 0)
allFR  	= cell2mat({r.allRes.dist2COM_corrFR}');
allFR  	= allFR - min(allFR, [], 2); % min = 0
allFR  	= allFR ./ range(allFR, 2); % max = 1

% define alphadata
alphaData                   = double(repmat(bDistCellResid, 1, size(allFR, 2))); % opaque
alphaData(alphaData == 0)   = 0.2; % transparent
alphaData(isnan(allFR))     = 0; % white

% sort tuning curves
[~, prefDistIdx]            = max(allFR, [], 2);
[~, sortIdx]                = sort(prefDistIdx, 'ascend');

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7, 16]);
axes('units', 'centimeters', 'position', [1, 1.75, 5, 12.75]);
imagesc(r.distnc.binCenters, 1:size(allFR, 1), allFR(sortIdx, :), ...
    'alphadata', alphaData(sortIdx, :));
cb = colorbar('southoutside');
cb.Units = 'centimeters';
cb.Position = [2.25, 14.75, 2.5, 0.3];
ylabel(cb, 'FR', ...
    'units', 'normalized', 'position', [0.5, 1.55]);
cb.Ticks = [0, 1];
cb.TickLabels = {'min', 'max'};
xl = xlabel('Distance (vu)');
yl = ylabel('Cell index', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set(gca, ...
    'xtick', [min(r.distnc.binCenters), max(r.distnc.binCenters)], ...
    'ytick', [min(sortIdx), max(sortIdx)], ...
    'ticklength', [0, 0]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
colormap jet;
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ELRPDCells_distanceTuning_20200922'), '-dtiff', '-r600');

%% identify distance cells based on 2D bearing-distance plot

% report
fprintf('\nEvaluation of bearing-distance tuning using the 2D bearing-distance plots.\n');
rng(r.param.myRNG);

% preallocate and loop through cells
allFiringRate               = cell(size(r.allRes, 1), 1); % all firing-rate maps
allFiringRateSurro          = cell(size(r.allRes, 1), 1);
allBearingDistanceField     = cell(size(r.allRes, 1), 1); % all bearing-distance fields
allMaxFieldSize             = nan(size(r.allRes, 1), 1); % size of all bearing-distance fields
allMaxFieldSizeSurro        = nan(size(r.allRes, 1), r.bearDist.numSurrogates);
for iCell = 1:size(r.allRes, 1)
    
    % report
    fprintf('Cell: %d (of %d). Index: %d, %d, %d.\n', iCell, size(r.allRes, 1), r.allRes(iCell).idx);
    
    %% empirical distance*bearing map
    
    % assemble relevant data
    bTmpMask            = r.allRes(iCell).bMask4Analysis & r.allRes(iCell).bAddMask4Analysis; % same data as for the other analysis
    thisFR              = r.allRes(iCell).FR(bTmpMask);
    thisDist2COM        = r.allRes(iCell).dist2COM(bTmpMask);
    thisEgoLocDir       = r.allBehLRP(iCell).egoLocDir(bTmpMask);
    
    % create 2D firing-rate map
    cfg                     = [];
    cfg.FR                  = thisFR;
    cfg.dist2COM            = thisDist2COM;
    cfg.egoLocDir           = thisEgoLocDir;
    cfg.distanceBinEdges    = r.bearDist.distanceBinEdges;
    cfg.distanceBinCenters  = r.bearDist.distanceBinCenters;
    cfg.bearingBinEdges     = r.bearDist.bearingBinEdges;
    cfg.bearingBinCenters   = r.bearDist.bearingBinCenters;    
    cfg.bSmooth             = true;
    cfg.smoothFac           = r.bearDist.smoothFac;
    cfg.smoothType          = r.bearDist.smoothType;
    empOut                  = LK_EstimateBearingDistanceMap_20210515(cfg);
    
    %% surrogate bearing-distance maps
    
    % random circular shifts of firing rates
    randShifts      = transpose(datasample(1:numel(thisFR) - 1, r.bearDist.numSurrogates, 'replace', false));
    
    % preallocate and loop through surrogates
    firingRateSurro = nan(size(empOut.firingRate, 1), size(empOut.firingRate, 2), r.bearDist.numSurrogates);
    for iSurro = 1:r.bearDist.numSurrogates
        
        % create surrogate firing rates
        surroThisFR             = circshift(thisFR, randShifts(iSurro));
        
        % create 2D distance-bearing ratemap
        cfg                     = [];
        cfg.FR                  = surroThisFR;
        cfg.dist2COM            = thisDist2COM;
        cfg.egoLocDir           = thisEgoLocDir;
        cfg.distanceBinEdges    = r.bearDist.distanceBinEdges;
        cfg.distanceBinCenters  = r.bearDist.distanceBinCenters;
        cfg.bearingBinEdges     = r.bearDist.bearingBinEdges;
        cfg.bearingBinCenters   = r.bearDist.bearingBinCenters;
        cfg.bSmooth             = true;
        cfg.smoothFac           = r.bearDist.smoothFac;
        cfg.smoothType          = r.bearDist.smoothType;
        surroOut                = LK_EstimateBearingDistanceMap_20210515(cfg);
        
        % collect surrogate values
        firingRateSurro(:, :, iSurro)   = surroOut.firingRate;
    end
    
    % collect results across cells
    allFiringRate{iCell, 1}    	= empOut.firingRate;
    allFiringRateSurro{iCell}  	= firingRateSurro;
    
    %% identify bearing*distance fields
    
    % empirical bearing-distance field
    cfgF                    = [];
    cfgF.firingRate         = allFiringRate{iCell};
    cfgF.firingRateSurro    = allFiringRateSurro{iCell};
    empField                = LK_DetectBearingDistanceFields_20201011(cfgF);
        
    % surrogate field strengths
    maxFieldSizeSurro   = nan(r.bearDist.numSurrogates, 1);
    surroIdx            = 1:r.bearDist.numSurrogates;
    for iSurro = 1:r.bearDist.numSurrogates
        
        % use hypothetical empirical and hypothetical surrogate data
        cfgF                    = [];
        cfgF.firingRate         = allFiringRateSurro{iCell}(:, :, surroIdx == iSurro);
        cfgF.firingRateSurro    = cat(3, allFiringRate{iCell}, allFiringRateSurro{iCell}(:, :, surroIdx ~= iSurro));
        surroField              = LK_DetectBearingDistanceFields_20201011(cfgF);
        
        % maximum surrogate field size
        if ~isempty(surroField.maxFieldSize)
            maxFieldSizeSurro(iSurro, 1)    = surroField.maxFieldSize;
        end
    end
    
    % collect across cells
    allBearingDistanceField{iCell, 1}   = empField.bearingDistanceField;
    allMaxFieldSize(iCell, 1)           = empField.maxFieldSize;
    allMaxFieldSizeSurro(iCell, :)      = maxFieldSizeSurro;
    
    %% figure
    
    % p-value for bearing-distance field
    pval4plot                   = 1 - sum(empField.maxFieldSize > maxFieldSizeSurro) / sum(~isnan(maxFieldSizeSurro));
    
    % create figure for bearing-distance maps
    bd                          = [];
    bd.visible                  = 'off';
    bd.bearingBinEdges          = r.bearDist.bearingBinEdges;
    bd.bearingBinCenters        = r.bearDist.bearingBinCenters;
    bd.distanceBinEdges         = r.bearDist.distanceBinEdges;
    bd.distanceBinCenters       = r.bearDist.distanceBinCenters;
    bd.firingRate               = allFiringRate{iCell};
    bd.bearingDistanceField     = allBearingDistanceField{iCell};
    bd.thisSpike                = r.allRes(iCell).thisSpike;
    bd.sr                       = r.allRes(iCell).t.par.sr;
    bd.nspk                     = sum(r.allRes(iCell).numSpikes(bTmpMask));
    bd.figTitle                 = {'Bearing-distance field', ...
        ['\itP\rm = ', num2str(pval4plot, '%.2f'), ' (', r.allRes(iCell).wireRegion, ')']};
    if pval4plot < 0.01
        bd.figTitle             = {'Bearing-distance field', ...
            ['\itP\rm < 0.01 (', r.allRes(iCell).wireRegion, ')']};
    elseif pval4plot < 0.05 && round(pval4plot, 2) == 0.05
        bd.figTitle             = {'Bearing-distance field', ...
            ['\itP\rm < 0.05 (', r.allRes(iCell).wireRegion, ')']};
    end
    figHandle   = LK_PlotBearingDistanceMap_20201011(bd);
    
    % save figure
    set(figHandle, 'PaperPositionMode', 'auto');
    print(figHandle, strcat(paths.save, subjects{r.allRes(iCell).idx(1)}, '_chan', num2str(r.allRes(iCell).idx(2)), ...
        '_Clus', num2str(r.allRes(iCell).idx(3)), '_Distance2BearingPlot_20201011'), '-dtiff', '-r300');
    close(figHandle);
end

%% distance-tuned ELRPD cells based on significant bearing-distance field

% distance-tuned ELRPD cells based on bearing-distance fields
maxFieldStrengthRanks   = sum(allMaxFieldSize > allMaxFieldSizeSurro, 2) ./ sum(~isnan(allMaxFieldSizeSurro), 2);
bDistCellField          = maxFieldStrengthRanks > 0.95;
fprintf('Number of distance-tuned ELRPD cells based on significant bearing-distance fields: %d.\n', sum(bDistCellField));

%% summary figure of bearing-distance fields

% concatenate all bearing-distance fields
catBearingDistanceField = nan(numel(r.bearDist.distanceBinCenters), numel(r.bearDist.bearingBinCenters), size(allBearingDistanceField, 1));
for iCell = 1:size(allBearingDistanceField, 1)
    catBearingDistanceField(:, :, iCell)    = allBearingDistanceField{iCell};
end

% different groups
groups  = {'allELRPDCells'};

% loop through groups
for iGroup = 1:numel(groups)
    
    % select data for this group of cells
    if strcmp(groups{iGroup}, 'allELRPDCells')
        sumBearingDistanceField     = sum(catBearingDistanceField(:, :, :), 3);
    end
    
    % summary 2D bearing-distance plot
    f = figure('units', 'centimeters', 'position', [2, 2, 7.5, 8], 'Color', [1, 1, 1]);
    
    % bearing-distance summary
    ax1 = axes('units', 'centimeters', 'position', [1.5, 1.75, 4, 4]);
    hold on;
    imagesc(r.bearDist.bearingBinCenters, r.bearDist.distanceBinCenters, sumBearingDistanceField, ...
        'alphadata', sumBearingDistanceField ~= 0);
    hold off;
    myCM = jet;
    colormap(myCM);
    cb = colorbar;
    ylabel(cb, 'Count', 'Rotation', 0, ...
        'units', 'normalized', 'position', [1.75, 0.5], ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    cb.Units = 'centimeters';
    cb.Position = [5.8, 6.05, 0.25, 1.25];
    cb.Limits = round([1, max(cb.Limits)], 1);
    cb.Ticks = round([1, max(cb.Limits)], 1);
    cb.Label.FontUnits = 'centimeters';
    cb.Label.FontSize = 0.4;
    xl = xlabel('Bearing');
    yl = ylabel('Distance (vu)', ...
        'units', 'normalized', 'position', [-0.05, 0.5]);
    idxValidY = find(sum(sumBearingDistanceField, 2) > 0, 1, 'first'):find(sum(sumBearingDistanceField, 2) > 0, 1, 'last') + 1;
    set(gca, ...
        'xlim', [min(r.bearDist.bearingBinEdges), max(r.bearDist.bearingBinEdges)], 'xtick', linspace(min(r.bearDist.bearingBinEdges), max(r.bearDist.bearingBinEdges), 5), ...
        'xticklabel', {'B', 'L', 'A', 'R', ''}, ...
        'ylim', [0, max(r.bearDist.distanceBinEdges(idxValidY))], 'ytick', [min(r.bearDist.distanceBinEdges), max(r.bearDist.distanceBinEdges(idxValidY))], ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    set([gca, xl, yl], ...
        'Fontunits', 'centimeters', 'Fontsize', 0.5);
    
    % distance summary
    ax2 = axes('units', 'centimeters', 'position', [5.8, 1.75, 1.25, 4]);
    sumDist = sum(sumBearingDistanceField, 2);
    plot(sumDist(idxValidY), r.bearDist.distanceBinCenters(idxValidY), '-', ...
        'Color', [0.5, 0.5, 0.5]);
    set(gca, ...
        'xlim', [0, max(sumDist(idxValidY))], 'xtick', [], ...
        'ylim', [0, max(r.bearDist.distanceBinEdges(idxValidY))], 'ytick', [], ...
        'box', 'off');
    
    % bearing summary
    ax3 = axes('units', 'centimeters', 'position', [1.5, 6.05, 4, 1.25]);
    sumBearing = sum(sumBearingDistanceField, 1);
    plot(r.bearDist.bearingBinCenters, sumBearing, '-', ...
        'Color', [0.5, 0.5, 0.5]);
    set(gca, ...
        'xlim', [min(r.bearDist.bearingBinEdges), max(r.bearDist.bearingBinEdges)], 'xtick', [], ...
        'ylim', [0, max(sumBearing)], 'ytick', [], ...
        'ytick', [], ...
        'box', 'off');
    
    % save figure
    set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    print(f, strcat(r.paths.save, 'SummaryDistance2BearingPlot_', groups{iGroup}, '_20201010'), '-dtiff', '-r300');
end

%% describe extent of bearing-distance fields

% xy-extent of bearing-distance fields
xyExtentFR          = nan(size(r.allRes, 1), 2); % bearing and distance extent of firing-rate maps
xyExtentField       = nan(size(r.allRes, 1), 2); % bearing and distance extent of bearing-distance fields
xyExtentFieldPerc   = nan(size(r.allRes, 1), 2);
extentFieldPerc     = nan(size(r.allRes, 1), 1); % total size of bearing-distance fields in percent
for iCell = 1:size(r.allRes, 1)
    
    % xy extent of entire firing-rate map
    xyExtentFR(iCell, :)        = [sum(any(~isnan(allFiringRate{iCell}), 1), 2), ...
        sum(any(~isnan(allFiringRate{iCell}), 2), 1)];
    
    % xy extent of bearing-distance field
    xyExtentField(iCell, :)     = [sum(any(allBearingDistanceField{iCell}, 1), 2), ...
        sum(any(allBearingDistanceField{iCell}, 2), 1)];
    
    % xy extent of bearing-distance field, in percent
    xyExtentFieldPerc(iCell, :) = xyExtentField(iCell, :) ./ xyExtentFR(iCell, :);
    
    % entire extent of bearing-distance field, in percent
    extentFieldPerc(iCell, 1)   = sum(allBearingDistanceField{iCell}(:)) / sum(~isnan(allFiringRate{iCell}(:)));
end

% create figure for bearing-distance extent of bearing-distance fields
f = figure('units', 'centimeters', 'position', [2, 2, 7, 7]);
hold on;
plot(100 .* xyExtentFieldPerc(:, 1), 100 .* xyExtentFieldPerc(:, 2), 'o', ...
    'Color', [0.7, 0.7, 0.7]);
plot(100 .* xyExtentFieldPerc(bDistCellField, 1), 100 .* xyExtentFieldPerc(bDistCellField, 2), 'o', ...
    'Color', [0.7, 0.7, 0.7], 'MarkerFaceColor', rgb('green'));
hold off;
xl = xlabel('Bearing extent (%)');
yl = ylabel('Distance extent (%)');
axis equal;
set(gca, ...
    'xlim', [0, 100], 'xtick', linspace(0, 100, 5), 'ylim', [0, 100], 'ytick', linspace(0, 100, 5), ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(r.paths.save, 'ExtentDistanceBearingFields_20201010'), '-dtiff', '-r300');

%% create figure for total size of bearing-distance fields

% extent bins
xEdges  = 100 .* linspace(0, myceil(max(extentFieldPerc), 2), 8);
fprintf('The total extent of significant bearing-distance fields ranges between %.3f and %.3f %%.\n', 100 * min(extentFieldPerc(bDistCellField)), 100 * max(extentFieldPerc(bDistCellField)));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 4, 4]);
axes('units', 'centimeters', 'position', [1, 1.6, 2.5, 2.1]);
hold on;
histogram(100 .* extentFieldPerc(~bDistCellField), xEdges, ...
    'FaceColor', [0.7, 0.7, 0.7], 'facealpha', 0.5);
histogram(100 .* extentFieldPerc(bDistCellField), xEdges, ...
    'FaceColor', rgb('green'), 'facealpha', 0.5);
xl = xlabel('Extent (%)');
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
tmpAx = get(gca);
set(gca, ...
    'xlim', [min(xEdges), max(xEdges)], 'xtick', [min(xEdges), max(xEdges)], ...
    'ytick', tmpAx.YLim, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(r.paths.save, 'ExtentDistanceBearingFieldsTotal_20201010'), '-dtiff', '-r300');
