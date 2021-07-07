%==========================================================================
% This script analyzes to what extent ELRPD-cell tuning relies on goal
% tuning.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));
parpool('local', 8);

% variables
param                   = [];
param.myRNG             = 9999;
param.clusterName       = 'Cluster4Analysis_HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768; % maximum absolute yaw value in the arena
param.newTimeRes        = 0.1; % temporal resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates   	= 1001; % number of surrogates
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
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sums of squares

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.save      = strcat('E:\OpenField\ELRPDCellGoalTuning_20201122\20201130_Dir', num2str(direc.angularRes), ...
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

%% load previous ELRPD cell results

% load previously saved ELRPD cell results
ELRPD  = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'locDir');

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByRefAna      = []; % excluded due to low number of spikes

% behavioral bookkeeping
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
    allBeh  = cat(1, allBeh, outBeh);
    
    %% process behavioral data with respect to local reference points
    
    % define goals (8 correct object locations; 8 average response
    % locations)
    objNames    = unique(outBeh.trials(:, 2));
    objLocs     = nan(size(objNames, 1), 2); % object locations (xy)
    respLocs    = nan(size(objNames, 1), 2); % response locations (xy)
    for iObj = 1:numel(objNames)
        logIdx              = outBeh.trials(:, 2) == objNames(iObj) & ~isnan(outBeh.trials(:, 9));
        objLocs(iObj, :)    = mean(outBeh.trials(logIdx, 9:10), 1);
        respLocs(iObj, :)   = mean(outBeh.trials(logIdx, 11:12), 1);
    end
    
    % use the object- and response-locations as LRPs (subject-specific)
    subjLRPs     = [objLocs; respLocs]; % object locations, response locations
    subjObjNames = [objNames; objNames]; % object identity
    
    % process behavior with respect to local reference points
    egoLocDirBins   = cell(size(subjLRPs, 1), 1);
    parfor iLRP = 1:size(subjLRPs, 1)
        
        % calculate egocentric directions
        alloLocDir  = atan2(subjLRPs(iLRP, 2) - outBeh.behinfo(:, 3), subjLRPs(iLRP, 1) - outBeh.behinfo(:, 2)); % atan2(y, x)
        egoLocDir   = angdiff(outBeh.yaws, alloLocDir); % angdiff(x, y) returns the angular difference delta = y - x
        
        % discretize egocentric directions
        egoLocDirBins{iLRP, 1}  = discretize(egoLocDir, locDir.angleEdges);
    end
    
    %% calculate behavioral data for surrogate reference points of interest
    
    % create random angles for shifting the goal locations
    randAngles  = rand(param.numSurrogates, 1) .* 2 .* pi;
    
    % create surrogate behavioral data
    subjLRPsSurro       = nan(size(subjLRPs, 1), size(subjLRPs, 2), param.numSurrogates); % surrogate LRPs; e.g., size = 16 x 2 x 101
    egoLocDirBinsSurro  = cell(size(subjLRPs, 1), param.numSurrogates); % surrogate egocentric directions; e.g., size = 16 x 101
    for iSurro = 1:param.numSurrogates
        
        % circularly shift the locations by a common angle
        [th, rh]        = cart2pol(objLocs(:, 1), objLocs(:, 2));
        [x, y]          = pol2cart(th + randAngles(iSurro), rh);
        surroObjLocs    = [x, y];
        [th, rh]        = cart2pol(respLocs(:, 1), respLocs(:, 2));
        [x, y]          = pol2cart(th + randAngles(iSurro), rh);
        surroRespLocs   = [x, y];
        
        % use the object- and response-locations as LRPs (subject-specific)
        subjLRPsSurro(:, :, iSurro) = [surroObjLocs; surroRespLocs]; % object locations, response locations
        
        % process behavior with respect to local reference points
        parfor iLRP = 1:size(subjLRPsSurro, 1)
            
            % calculate egocentric directions
            alloLocDir  = atan2(subjLRPsSurro(iLRP, 2, iSurro) - outBeh.behinfo(:, 3), subjLRPsSurro(iLRP, 1, iSurro) - outBeh.behinfo(:, 2)); % atan2(y, x)
            egoLocDir   = angdiff(outBeh.yaws, alloLocDir); % angdiff(x, y) returns the angular difference delta = y - x
            
            % discretize egocentric directions
            egoLocDirBinsSurro{iLRP, iSurro}    = discretize(egoLocDir, locDir.angleEdges);
        end
    end
    
    %% mask for analysis timepoints (timepoints to be included in the analysis)
    
    % boolean that encodes whether specific time bins shall be included in
    % the analysis
    bMask4Analysis  = outBeh.bNavi | outBeh.bRotation;
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
        
        %% brain region of this wire
        
        % subject information
        s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion  = s.subjectdata.micro2macro{logIdx, 3};
        
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
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, param.clusterName, '.mat'));
        
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
                        
            % continue if this cluster was not included in the ELRPD cell
            % analysis or if it is not an ELRPD cell
            refIdx  = all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2);
            if sum(refIdx) ~= 1 || ~(ELRPD.allRes(refIdx).locdir_clusRank > 0.95)
                fprintf('\t\t- This cluster was not included in the ELRPD cell analysis or is not an ELRPD cell.\n');
                exByRefAna  = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
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
            
            %% ANOVA for modulation of firing rate by egocentric direction towards local reference points
            
            % for each local reference point, estimate modulation of firing
            % rate by egocentric direction towards this LRP
            resFile = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_locdir_emp.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp             = load(resFile);
                locdir_F        = tmp.locdir_F; % F value
                locdir_corrFR   = tmp.locdir_corrFR; % tuning curve
            else
                
                % preallocate
                locdir_F        = nan(size(subjLRPs, 1), 1); % e.g., size = 16 x 1
                locdir_corrFR   = nan(size(subjLRPs, 1), numel(locDir.angleCenters)); % e.g., size = 16 x 12
                for iLRP = 1:size(subjLRPs, 1)
                    
                    % compute ANOVAN
                    dt              = [];
                    dt.FR           = FR(bMask4Analysis);
                    dt.group        = [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis), egoLocDirBins{iLRP}(bMask4Analysis)];
                    dt.groupNames   = {'D'; 'P'; 'ELRPD'}; % predictors
                    dt.ANOVA        = ANOVA;
                    dt.bSkipA1      = 'yes';
                    empANOVA        = LK_CompANOVA_070719(dt);
                    
                    % output F value
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
            
            % create surrogates or load previously saved surrogates
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                 = load(resFile);
                locdir_F_surro     	= tmp.locdir_F_surro;
                locdir_corrFRsurro  = tmp.locdir_corrFRsurro;
            else
                
                % preallocate
                locdir_F_surro  	= nan(size(subjLRPsSurro, 1), param.numSurrogates);
                locdir_corrFRsurro  = nan(size(subjLRPsSurro, 1), numel(locDir.angleCenters), param.numSurrogates);
                
                % loop through surrogates
                parfor iSurro = 1:param.numSurrogates
                                                            
                    %% surrogate ANOVA for modulation of firing rate by egocentric direction towards local reference points
                    
                    % for each local reference point, estimate modulation
                    % of firing rate by egocentric direction towards this
                    % LRP
                    tmp_locdir_F_surro      = nan(size(subjLRPsSurro, 1), 1);
                    tmp_locdir_corrFRsurro  = nan(size(subjLRPsSurro, 1), numel(locDir.angleCenters));
                    for iLRP = 1:size(subjLRPsSurro, 1)
                        
                        % compute ANOVAN
                        dt              = [];
                        dt.FR           = FR(bMask4Analysis); % cave: use empirical firing rates
                        dt.group        = [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis), egoLocDirBinsSurro{iLRP, iSurro}(bMask4Analysis)];
                        dt.groupNames   = {'D'; 'P'; 'ELRPD'}; % predictors
                        dt.ANOVA        = ANOVA;
                        dt.bSkipA1      = 'yes';
                        surroANOVA      = LK_CompANOVA_070719(dt);
                        
                        % get output F value
                        tmp_locdir_F_surro(iLRP, 1)         = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'ELRPD'), strcmp(surroANOVA.tbl(1, :), 'F')));
                        
                        % tuning curve for ELRPD
                        tmpFR                               = surroANOVA.m{strcmp(surroANOVA.groupNames, 'ELRPD')}(:, 1);
                        idx                                 = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'ELRPD')});
                        tmp_locdir_corrFRsurro(iLRP, idx) 	= tmpFR;
                    end
                    
                    % collect across surrogate rounds
                    locdir_F_surro(:, iSurro)               = tmp_locdir_F_surro;
                    locdir_corrFRsurro(:, :, iSurro)        = tmp_locdir_corrFRsurro;
                end
                
                % save surrogates
                save(resFile, ...
                    'locdir_F_surro', 'locdir_corrFRsurro');
            end
            
            %% statistical evaluation of egocentric bearing tuning
            
            % rank of ELRPD F-value
            locdir_Frank    = sum(locdir_F > locdir_F_surro, 2) ./ sum(~isnan(locdir_F_surro), 2);
                        
            %% collect information for this unit
            
            % basics
            unitRes                             = [];
            unitRes.idx                         = [iSub, iWire, iClus];
            unitRes.numSpikes                   = numSpikes;
            unitRes.bMask4Analysis              = bMask4Analysis;
            unitRes.wireRegion                  = wireRegion;
            
            % ELRPD
            unitRes.locdir_F                    = locdir_F;
            unitRes.locdir_F_surro              = locdir_F_surro;
            unitRes.locdir_Frank                = locdir_Frank;
            unitRes.locdir_corrFR               = locdir_corrFR;
            unitRes.locdir_corrFRsurro          = locdir_corrFRsurro;
                        
            % save unit result
            save(strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_unitRes'), 'unitRes');
            
            % collect across units  
            allRes  = cat(1, allRes, unitRes);
            
            %% figure
            
            % estimate smallest p-value for any of the goal tunings
            minP                        = 1 - max(locdir_Frank);
            
            % data for figure
            dt                          = [];
            dt.visible                  = 'off';
            dt.locDir                   = locDir;
            dt.locDir.LRPs              = subjLRPs;
            dt.locdir_Frank             = locdir_Frank;
            dt.locdir_corrFR            = locdir_corrFR;
            dt.thisSpike                = thisSpike; % spike waveforms
            dt.t                        = t; % includes sampling rate
            dt.numSpikes                = numSpikes;
            dt.bMask4Analysis           = bMask4Analysis;
            dt.figTitle                 = {'Egocentric goal tuning', ...
                ['\itP\rm_{min} = ', num2str(minP, '%.2f'), ' (', wireRegion, ')']};
            if minP < 0.05 && round(minP, 2) == 0.05
                dt.figTitle   	= {'Egocentric goal tuning', ...
                    ['\itP\rm_{min} < 0.05 (', wireRegion, ')']};
            elseif minP < 0.01
                dt.figTitle  	= {'Egocentric goal tuning', ...
                    ['\itP\rm_{min} < 0.01 (', wireRegion, ')']};
            end
            
            % plot figure
            f = LK_OF_PlotGoalTuning_20201122(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_GoalTuning'), '-dtiff', '-r450');
            
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

%% evaluation of goal tuning

% empirical F and surrogate F values for each cell
allF        = cell2mat({r.allRes.locdir_F})';
allFsurro   = nan(size(r.allRes, 1), size(r.allRes(1).locdir_F_surro, 1), size(r.allRes(1).locdir_F_surro, 2));
for iCell = 1:size(r.allRes, 1)
    allFsurro(iCell, :, :)  = r.allRes(iCell).locdir_F_surro;
end

% select specific type of objects to examine
selObj      = [true(8, 1); false(8, 1)] == true; % true ==> object locations; false ==> response locations

% choose maximum across objects
maxF        = max(allF(:, selObj), [], 2);
maxFsurro   = squeeze(max(allFsurro(:, selObj, :), [], 2));

% test how often the empirical maximum is above the surrogate maxima
maxFrank    = sum(maxF > maxFsurro, 2) ./ sum(~isnan(maxFsurro), 2);

% perform t-test of F ranks against chance level (0.5)
[~, p, ~, stats]    = ttest(maxFrank, 0.5);
fprintf('T-test for "maxFrank"-values against 0.5 chance: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
