%==========================================================================
% This script analyzes egocentric bearing cells (= ELRPD cells) using an
% ANOVA with a reduced number of predictors.
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
param.maxAbsYaw         = 32768; % maximum absolute yaw value in the arena
param.newTimeRes        = 0.1; % new time resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates   	= 101; % number of surrogates
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
%- for location-specific directional tuning
locDir.areaAroundLRP    = 3333; % 1/3 of arena diameter
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per behavioral bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sums of squares

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.arena     = 'E:\OpenField\Arena\';
paths.save      = strcat('E:\OpenField\ELRPDCellAnalysis_20200623\20200617', ...
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

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByVisInsp     = []; % excluded by visual inspection
exByNumSpikes   = []; % excluded due to low number of spikes

% for behavioral control analyses
allBeh          = [];
allBehLRP       = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(param.myRNG);
    
    % subject-specific save path
    subjSavePath    = strcat(paths.save, subjects{iSub}, '\');
    if ~exist(subjSavePath, 'dir')
        mkdir(subjSavePath); % create folder
    end
    
    %% process behavioral data
    
    % process behavior
    behFile     = strcat(subjSavePath, 'outBeh.mat');
    if exist(behFile, 'file')
        
        % load previously processed behavioral data
        fprintf('Loading behavioral data from "%s".\n', behFile);
        tmp                 = load(behFile);
        outBeh              = tmp.outBeh;
    else
        
        % process behavior
        dt                  = [];
        dt.paths            = paths;
        dt.subject          = subjects{iSub};
        dt.direc            = direc;
        dt.place            = place;
        dt.maxAbsYaw        = param.maxAbsYaw;
        dt.maxNumTrials     = param.maxNumTrials;        
        dt.newTimeRes       = param.newTimeRes;
        dt.naviCutoff       = param.naviCutoff;
        dt.naviSmooth       = param.naviSmooth;
        dt.rotSmooth        = param.rotSmooth;
        outBeh              = LK_EvalBehOverall_20200623(dt);
        
        % use new timeline for the analyses
        outBeh              = outBeh.new;
    
        % save output
        save(behFile, 'outBeh');
    end
    
    % collect behavioral output for all subjects
    allBeh              = cat(1, allBeh, outBeh);
    
    %% process behavioral data for all local reference points
        
    % process behavior with respect to local reference points
    behLRPFile  = strcat(subjSavePath, 'outBehLRP.mat');
    if exist(behLRPFile, 'file')
        
        % load previously processed behavioral data
        tmp      	= load(behLRPFile);
        outBehLRP  	= tmp.outBehLRP;
    else
        outBehLRP  	= cell(size(locDir.LRPs, 1), 1);
        parfor iLRP = 1:size(locDir.LRPs, 1)
            
            % process behavior with specific local reference point
            dt                  = [];
            dt.paths            = paths;
            dt.subject          = subjects{iSub};
            dt.direc            = direc;
            dt.place            = place;
            dt.locDir           = locDir;
            dt.locDir.refPoint  = locDir.LRPs(iLRP, :); % specific LRP
            dt.maxAbsYaw        = param.maxAbsYaw;
            dt.maxNumTrials     = param.maxNumTrials;
            dt.newTimeRes       = param.newTimeRes;
            dt.naviCutoff       = param.naviCutoff;
            dt.naviSmooth       = param.naviSmooth;
            dt.rotSmooth        = param.rotSmooth;
            tmpOut              = LK_EvalBehOverall_20200623(dt);
            
            % use new timeline for the analyses
            outBehLRP{iLRP, 1}  = tmpOut.new;
        end
        
        % save output
        save(behLRPFile, 'outBehLRP');
    end
    
    % collect behavioral output for all subjects
    allBehLRP           = cat(1, allBehLRP, outBehLRP');
    
    %% mask for analysis timepoints
    
    % boolean that encodes whether specific time bins shall be included in
    % the analysis
    bMask4Analysis      = outBeh.bNavi | outBeh.bRotation;
    fprintf('Information: %.3f%% of the "active time" is included into the analysis.\n', ...
        sum(outBeh.durations(bMask4Analysis)) / sum(outBeh.durations(outBeh.bActive)) * 100);
    
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
        c4a     = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', param.clusterName, '.mat'));
        
        % load behavioral times of spikes
        ccb     = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('Size of "cluster_class" not congruent with "cluster_class_behtime');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
                        
            % continue if you decided that this cluster is insufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.Cluster}) == iClus).Decision, 'no')
                fprintf('\t\t- You decided not to analyse this cluster.\n');
                exByVisInsp     = cat(1, exByVisInsp, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% number of spikes and FR per time bin
            
            % number of spikes per timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time-bin
            
            % FR per timepoint
            FR  = numSpikes ./ outBeh.durations;
            
            % if the general firing rate is too low (<0.1 Hz), continue
            if (sum(numSpikes(bMask4Analysis)) / sum(outBeh.durations(bMask4Analysis))) < 0.1
                fprintf('\t\t- The mean firing rate of this cluster is less than 0.1 Hz, thus skipping.\n');
                exByNumSpikes   = cat(1, exByNumSpikes, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
                        
            %% ANOVA for modulation of firing rate by egocentric direction towards local reference points
            
            % for each local reference point (LRP), estimate modulation of
            % firing rate by egocentric direction towards this point
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_locdir_emp.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp             = load(resFile);
                locdir_F        = tmp.locdir_F; % F-value
                locdir_corrFR   = tmp.locdir_corrFR; % tuning curve
            else
                
                % preallocate
                locdir_F        = nan(size(locDir.LRPs, 1), 1);
                locdir_corrFR   = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters));
                for iLRP = 1:size(locDir.LRPs, 1)
                    
                    % compute ANOVAN
                    dt              = [];
                    dt.FR           = FR(bMask4Analysis);
                    dt.group        = [outBeh.xyBin(bMask4Analysis), outBehLRP{iLRP}.egoLocDirBin(bMask4Analysis)];
                    dt.groupNames   = {'P'; 'ELRPD'}; % predictors
                    dt.ANOVA        = ANOVA;
                    dt.bSkipA1      = 'yes';
                    empANOVA        = LK_CompANOVA_070719(dt);
                    
                    % get output F values
                    locdir_F(iLRP, 1)           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'ELRPD'), strcmp(empANOVA.tbl(1, :), 'F')));
                    
                    % tuning curve for ELRPD
                    tmpFR                       = empANOVA.m{strcmp(empANOVA.groupNames, 'ELRPD')}(:, 1);
                    idx                         = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'ELRPD')});
                    locdir_corrFR(iLRP, idx)    = tmpFR;
                end
                
                % save result
                save(strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_locdir_emp'), ...
                    'locdir_F', 'locdir_corrFR');
            end
            
            %% surrogates
            
            % random numbers for circular shift of firing rates
            randShifts  = transpose(datasample(1:numel(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
            
            % create surrogates or load previously saved surrogates
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                     = load(resFile);
                locdir_F_surro          = tmp.locdir_F_surro; % surrogate F values
                all_locdir_corrFRsurro  = tmp.all_locdir_corrFRsurro; % surrogate tuning curves
            else
                
                % preallocate
                locdir_F_surro          = nan(size(locDir.LRPs, 1), param.numSurrogates);
                all_locdir_corrFRsurro  = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters), param.numSurrogates);
                
                % loop through surrogates
                parfor iSurro = 1:param.numSurrogates
                    
                    % create surrogate firing rate by circularly shifting
                    % firing rate relative to the behavioral data
                    FRsurro         = circshift(FR(bMask4Analysis), randShifts(iSurro, 1)); % mask by analysis time window
                                        
                    %% surrogate ANOVA for modulation of firing rate by egocentric direction towards local reference points
                    
                    % for each local reference point, estimate modulation
                    % of firing rate by egocentric direction towards this
                    % LRP
                    tmp_locdir_F_surro  = nan(size(locDir.LRPs, 1), 1);
                    locdir_corrFRsurro  = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters));
                    for iLRP = 1:size(locDir.LRPs, 1)
                        
                        % compute ANOVAN
                        dt              = [];
                        dt.FR           = FRsurro; % you already masked the FR by the analysis period (see above)
                        dt.group        = [outBeh.xyBin(bMask4Analysis), outBehLRP{iLRP}.egoLocDirBin(bMask4Analysis)];
                        dt.groupNames   = {'P'; 'ELRPD'};
                        dt.ANOVA        = ANOVA;
                        dt.bSkipA1      = 'yes';
                        surroANOVA      = LK_CompANOVA_070719(dt);
                        
                        % get output F value
                        tmp_locdir_F_surro(iLRP, 1)     = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'ELRPD'), strcmp(surroANOVA.tbl(1, :), 'F')));
                        
                        % tuning curve for ELRPD
                        tmpFR                           = surroANOVA.m{strcmp(surroANOVA.groupNames, 'ELRPD')}(:, 1);
                        idx                             = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'ELRPD')});
                        locdir_corrFRsurro(iLRP, idx) 	= tmpFR;
                    end
                    
                    % collect across surrogate rounds
                    locdir_F_surro(:, iSurro)               = tmp_locdir_F_surro;
                    all_locdir_corrFRsurro(:, :, iSurro)    = locdir_corrFRsurro;                    
                end
                
                % save surrogates
                save(resFile, ...
                    'locdir_F_surro', 'all_locdir_corrFRsurro');
            end
            
            %% allocentric direction tuning at each candidate LRP (+ surrounding area)
            
            % perform analysis or load previously saved results
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_dir_corrFR_perLRP.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                 = load(resFile);
                dir_F_perLRP        = tmp.dir_F_perLRP;
                dir_corrFR_perLRP   = tmp.dir_corrFR_perLRP;
            else
                
                % loop through LRPs
                dir_F_perLRP        = nan(size(locDir.LRPs, 1), 1);
                dir_corrFR_perLRP   = nan(size(locDir.LRPs, 1), size(direc.angleCenters, 1));
                for iLRP = 1:size(locDir.LRPs, 1)
                    
                    % select time points when subject was in the vicinity
                    % of this LRP
                    bMask4AnalysisLRP   = bMask4Analysis & ...
                        outBeh.behinfo(:, 2) >= (locDir.LRPs(iLRP, 1) - locDir.areaAroundLRP) & ...
                        outBeh.behinfo(:, 2) <= (locDir.LRPs(iLRP, 1) + locDir.areaAroundLRP) & ...
                        outBeh.behinfo(:, 3) >= (locDir.LRPs(iLRP, 2) - locDir.areaAroundLRP) & ...
                        outBeh.behinfo(:, 3) <= (locDir.LRPs(iLRP, 2) + locDir.areaAroundLRP);
                    
                    % compute ANOVAN for this data part
                    try
                        dt              = [];
                        dt.FR           = FR(bMask4AnalysisLRP);
                        dt.group        = [outBeh.yawBin(bMask4AnalysisLRP), outBeh.xyBin(bMask4AnalysisLRP)];
                        dt.groupNames   = {'D'; 'P'};
                        dt.ANOVA        = ANOVA;
                        dt.bSkipA1      = 'yes';
                        empANOVA        = LK_CompANOVA_070719(dt);
                        
                        % get output F value
                        dir_F_perLRP(iLRP, 1)           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'D'), strcmp(empANOVA.tbl(1, :), 'F')));
                        
                        % corrected firing rate for direction
                        tmpFR                           = empANOVA.m{strcmp(empANOVA.groupNames, 'D')}(:, 1);
                        idx                             = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')});
                        dir_corrFR_perLRP(iLRP, idx)    = tmpFR;
                    catch
                        fprintf('\t\t--- Location-specific analysis of direction tuning was not possible for LRP (%d/%d).\n', ...
                            locDir.LRPs(iLRP, 1), locDir.LRPs(iLRP, 2));
                    end
                end
                
                % save output
                save(resFile, 'dir_corrFR_perLRP', 'dir_F_perLRP');
            end
            
            %% statistical evaluation of ELRPD tuning
            
            % rank of ELRPD F-value
            locdir_Frank            = sum(locdir_F > locdir_F_surro, 2) ./ sum(~isnan(locdir_F_surro), 2); % e.g., 144 x 1
            
            % test whether the modulation map is significant in its total
            % (cluster-based permutation test)
            dt                      = [];
            dt.locDir               = locDir;
            dt.locRefPoints         = locDir.LRPs;
            dt.locRefPoints2Exclude = locDir.LRPs2Exclude;
            dt.locdir_F             = locdir_F;
            dt.locdir_F_surro       = locdir_F_surro;
            dt.locdir_Frank         = locdir_Frank;
            dt.locdir_corrFR        = locdir_corrFR;
            dt.numSurrogates        = param.numSurrogates;
            outCBPT                 = LK_CBPT_locdir_20200623(dt);
            
            %% collect information for this unit
            
            % basics
            unitRes                             = [];
            unitRes.idx                         = [iSub, iWire, iClus];
            unitRes.bMask4Analysis              = bMask4Analysis;
            unitRes.wireRegion                  = wireRegion;
            unitRes.wireMNI                     = wireMNI;
            
            % ELRPD
            unitRes.locdir_F                    = locdir_F;
            unitRes.locdir_F_surro              = locdir_F_surro;
            unitRes.locdir_Frank                = locdir_Frank;
            unitRes.locdir_corrFR               = locdir_corrFR;
            unitRes.locdir_clusRank             = outCBPT.locdir_clusRank;
            unitRes.locdir_maxSumF_emp          = outCBPT.maxSumF_emp;
            unitRes.locdir_maxSumF_surro        = outCBPT.maxSumF_surro;
            unitRes.locdir_largestCluster       = outCBPT.largestCluster;
            unitRes.locdir_COMxy                = outCBPT.COMxy;
            unitRes.locdir_COMxyClosestLRP      = outCBPT.COMxyClosestLRP;
            
            % direction, separately for each LRP
            unitRes.dir_F_perLRP                = dir_F_perLRP;
            unitRes.dir_corrFR_perLRP           = dir_corrFR_perLRP;
            
            % save unit-result
            save(strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_unitRes'), 'unitRes');
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);
            
            %% figure: ELRPD tuning
            
            % data for figure
            dt                  	= [];
            dt.visible              = 'off';
            dt.locDir           	= locDir;
            dt.outCBPT           	= outCBPT;
            dt.locdir_corrFR     	= locdir_corrFR;
            dt.direc             	= direc;
            dt.dir_corrFR_perLRP  	= dir_corrFR_perLRP;
            dt.thisSpike        	= thisSpike; % spike waveforms
            dt.t.par.sr         	= t.par.sr; % sampling rate
            dt.nspk                 = sum(numSpikes(bMask4Analysis)); % number of spikes
            dt.figTitle         	= {'Reference field and point', ...
                ['\itP\rm = ', num2str(1 - outCBPT.locdir_clusRank, '%.2f'), ' (', wireRegion, ')']};
            if (1 - outCBPT.locdir_clusRank) < 0.05 && round(1 - outCBPT.locdir_clusRank, 2) == 0.05
                dt.figTitle   	= {'Reference field and point', ...
                    ['\itP\rm < 0.05 (', wireRegion, ')']};
            elseif (1 - outCBPT.locdir_clusRank) < 0.01
                dt.figTitle  	= {'Reference field and point', ...
                    ['\itP\rm < 0.01 (', wireRegion, ')']};
            end
            
            % plot figure
            f = LK_OF_PlotELRPDModulation_20210408(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_ELRPD'), '-dtiff', '-r450');
            close(f);
            
            %% close all open figures
            close all;
            toc            
        end
    end
end

%% save all results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'), ...
    'allRes', 'allBeh', 'allBehLRP', ...
    'param', 'direc', 'locDir', 'paths', 'subjects');

% all cells' indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');
fprintf('Total number of cells: %d.\n', size(allUnitIdx, 1));

%% ELRPD cells

% identify ELRPD cells and report
all_locdir_clusRank = cell2mat({r.allRes.locdir_clusRank}');
bELRPDCell          = all_locdir_clusRank > 0.95; % >0.95 corresponds to P<0.05
fprintf('\nNumber of ELRPD cells: %d (P = %.3f).\n', sum(bELRPDCell), ...
    myBinomTest(sum(bELRPDCell), numel(bELRPDCell), 0.05));

%% compare to ELRPD cell analysis using the 3way ANOVA framework

% report
fprintf('Comparing this result to the ELRPD-cell result obtained from the 3way ANOVA:\n');

% load previous results
res3Way                     = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', 'allRes');
all_locdir_clusRank_3Way    = cell2mat({res3Way.allRes.locdir_clusRank}');
bELRPDCell_3Way             = all_locdir_clusRank_3Way > 0.95;

% correlation
[rho, pval]     = corr(all_locdir_clusRank, all_locdir_clusRank_3Way, 'type', 'spearman');
fprintf('Spearman correlation between 2-way and 3-way ANOVA results: rho = %.3f, p = %.3f.\n', rho, pval);

% chi-squared test
n   = [sum(bELRPDCell & bELRPDCell_3Way), sum(bELRPDCell & ~bELRPDCell_3Way); ...
    sum(~bELRPDCell & bELRPDCell_3Way), sum(~bELRPDCell & ~bELRPDCell_3Way)];
[X, p]  = myChiSquareTest(n(1, 1), n(1, 2), n(2, 1), n(2, 2));
fprintf('Chi-squared test: X = %.3f, p = %.3f.\n', X, p);
