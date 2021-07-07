%==========================================================================
% This script analyzes direction cells and place-like cells.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath('E:\TreasureHunt\Functions\');
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.myRNG             = 777;
param.clusterName       = 'Cluster4Analysis_LK_20200416_113912';
param.newTimeRes        = 0.1; % temporal resolution (sec)
param.naviCutoff        = 0.1; % speed cutoff
param.naviSmooth        = 2;
param.rotationCutoff    = 0.01; % angular speed cutoff
param.rotationSmooth    = 2;
param.arenaXLim         = [420, 320];
param.arenaZLim         = [410, 310];
param.arenaCtr          = [mean(param.arenaXLim), mean(param.arenaZLim)]; % arena center (x/z)
param.arenaRadius       = 50; % arena radius
param.numSurrogates     = 101; % number of surrogates
%- for direction
direc                   = [];
direc.angularRes        = 30; % angular resolution
direc.angleEdges        = deg2rad(-180:direc.angularRes:180);
direc.angleCenters      = movmean(direc.angleEdges, 2, 'endpoints', 'discard');
%- for location
place                   = [];
place.locRes            = 6; % spatial resolution
place.idxTemplate       = flipud(reshape(1:place.locRes^2, place.locRes, place.locRes)); % template for location indexing
place.xEdges            = linspace(430, 310, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.zEdges            = linspace(420, 300, place.locRes + 1);
place.zCenters          = movmean(place.zEdges, 2, 'endpoints', 'discard');
place.cutoffSigBin      = 0.95;
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per behavioral bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sums of squares

% paths
paths     	= [];
paths.info  = 'E:\TreasureHunt\SubjectInformation_170818\';
paths.spike = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.beh   = 'E:\TreasureHunt\Beh_210420\';
paths.save  = strcat('E:\TreasureHunt\PlaceDirAnalysis_220520\20200623_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.zCenters)), '\', ...
    'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotationIfAbove', num2str(param.rotationCutoff), ...
    '_Smooth', num2str(param.naviSmooth), '\');
mkdir(paths.save);

% subjects
subjects    = {...
    'TH_001'; ...
    'TH_002'; ...
    'TH_003'; ...
    'TH_004'; ...
    'TH_005'; ...
    'TH_006'; ...
    'TH_007'; ...
    'TH_008'; ...
    'TH_009'; ...
    'TH_010'; ...
    'TH_011'; ...
    'TH_012'; ...
    };

%% save settings
save(strcat(paths.save, 'settings'));

%% preallocations

% main results for all cells
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByVisInsp     = []; % excluded by visual inspection
exByNumSpikes   = []; % excluded due to low number of spikes

% all behavioral data
allBeh          = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % set rng
    rng(param.myRNG);
    
    % available sessions
    sessions    = dir(fullfile(paths.spike, subjects{iSub}, 'session_*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % report
        fprintf('\n\nSUBJECT: %s, SESSION: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % subject-specific save path
        subjSavePath    = fullfile(paths.save, subjects{iSub}, sessions(iSess).name, '\');
        if ~exist(subjSavePath, 'dir')
            mkdir(subjSavePath); % create folder
        end
        
        %% behavioral data
        
        % process behavior
        behFile     = strcat(subjSavePath, 'outBeh.mat');
        if exist(behFile, 'file')
            
            % load previously saved behavioral data
            fprintf('Loading behavioral data from "%s".\n', behFile);
            tmp     = load(behFile);
            outBeh  = tmp.outBeh;
        else            
            % load preprocessed behavioral data
            tmp     = load(fullfile(paths.beh, subjects{iSub}, sessions(iSess).name, ['behInfo', num2str(1 / param.newTimeRes), 'Hz.mat']));
            outBeh  = tmp.behInfo.new; % use new timeline
            
            % discretize yaw values
            outBeh.yawBins  = discretize(outBeh.yaws, direc.angleEdges);
            
            % discretize xz locations
            outBeh.xBins    = numel(place.xEdges) - discretize(outBeh.xyz(:, 1), fliplr(place.xEdges)); % high values get low bin indices
            outBeh.zBins    = numel(place.zEdges) - discretize(outBeh.xyz(:, 3), fliplr(place.zEdges)); % high values get low bin indices
            outBeh.xzBins   = (outBeh.xBins - 1) .* numel(place.zCenters) + outBeh.zBins; % lower left corner to upper right corner
                        
            % save output
            save(behFile, 'outBeh');
        end
        
        % collect behavioral output for all subjects
        allBeh  = cat(1, allBeh, outBeh);       
        
        %% mask for analysis timepoints (timepoints to be included in the analysis)
        
        % boolean vector for navigation
        outBeh.bNavi        = outBeh.speed > param.naviCutoff;
        if isfield(param, 'naviSmooth') && param.naviSmooth > 0
            convVec     	= ones(param.naviSmooth / param.newTimeRes, 1); % convolution vector
            outBeh.bNavi    = conv(outBeh.bNavi, convVec, 'same') > 0;
        end
        
        % boolean vector for rotation
        outBeh.bRotation        = outBeh.angSpeed > param.rotationCutoff;
        if isfield(param, 'rotationSmooth') && param.rotationSmooth > 0
            convVec             = ones(param.rotationSmooth / param.newTimeRes, 1); % convolution vector
            outBeh.bRotation    = conv(outBeh.bRotation, convVec, 'same') > 0;
        end
        
        % time points to be included into the analysis
        bMask4Analysis   	= strcmp(outBeh.behInfo(:, 10), 'navigation') & (outBeh.bNavi | outBeh.bRotation);
        fprintf('Information: %.3f%% (%.3f sec) of the entire logfile time is included into the analysis.\n', ...
            sum(outBeh.durations(bMask4Analysis)) / sum(outBeh.durations) * 100, sum(outBeh.durations(bMask4Analysis)));
        
        %% wires to investigate
        
        % available microwires
        wires   = dir(fullfile(paths.spike, subjects{iSub}, sessions(iSess).name, 'chan*'));
        tmp     = split(transpose({wires.name}), 'n');
        [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
        wires   = wires(I);
        
        %% loop through wires
        for iWire = 1:size(wires, 1)
            
            % report
            fprintf('\tWire: %s.\n', wires(iWire).name);
            
            %% brain region of this wire
            
            % subject information
            s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata.mat'));
            logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
            wireRegion  = s.subjectdata.micro2macro{logIdx, 3}; % brain region
            wireMNI     = s.subjectdata.micro2macro{logIdx, 5}; % MNI coordinates
            
            %% spike times in behavioral time
            
            % load wave-clus output
            wirePath    = strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep);
            try
                t = load(strcat(wirePath, 'times_datacut.mat'));
            catch
                fprintf('\t- No wave-clus for this wire.\n');
                exByWC  = cat(1, exByWC, [iSub, iSess, iWire]); % bookkeeping
                continue;
            end
            
            % load decision whether to use clusters (based on inspection)
            c4a = load(strcat(wirePath, param.clusterName, '.mat'));
            
            % load behavioral times of spikes
            ccb = load(strcat(wirePath, 'cluster_class_behtime.mat'));
            
            % sanity check
            if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
                error('Size of "cluster_class" not congruent with "cluster_class_behtime".');
            end
            
            %% loop through clusters
            for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
                
                % report
                fprintf('\t\tCluster: %d.\n', iClus);
                tic
                
                % continue if you decided that this cluster is insufficient
                if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.Cluster}) == iClus).Decision, 'no')
                    fprintf('\t\t- You decided not to analyse this cluster.\n');
                    exByVisInsp = cat(1, exByVisInsp, [iSub, iSess, iWire, iClus]); % bookkeeping
                    continue;
                end
                
                % data for this cluster
                thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
                thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
                
                %% number of spikes and firing rate per time bin
                
                % number of spikes per timepoint
                numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time-bin
                
                % firing rate per timepoint
                FR  = numSpikes ./ outBeh.durations;
                
                % if the general firing rate is too low (<0.1 Hz), continue
                if (sum(numSpikes(bMask4Analysis)) / sum(outBeh.durations(bMask4Analysis))) < 0.1
                    fprintf('\t\t- The mean firing rate of this cluster is less than 0.1 Hz, thus skipping.\n');
                    exByNumSpikes   = cat(1, exByNumSpikes, [iSub, iSess, iWire, iClus]); % bookkeeping
                    continue;
                end
                                
                %% ANOVA: direction and place on firing rate
                
                % compute ANOVAN
                dt              = [];
                dt.FR           = FR(bMask4Analysis); % mask with analysis period
                dt.group        = [outBeh.yawBins(bMask4Analysis), outBeh.xzBins(bMask4Analysis)];
                dt.groupNames   = {'D'; 'P'}; % predictors
                dt.ANOVA        = ANOVA;
                empANOVA        = LK_TH_CompANOVA_210520(dt);
                
                % output F values
                dir_F           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'D'), strcmp(empANOVA.tbl(1, :), 'F')));
                loc_F           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'P'), strcmp(empANOVA.tbl(1, :), 'F')));
                
                %% tuning curves
                                
                % tuning curve for direction
                dir_corrFR      = nan(numel(direc.angleCenters), 1);
                tmpFR           = empANOVA.m{strcmp(empANOVA.groupNames, 'D')}(:, 1);
                idx             = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')});
                dir_corrFR(idx) = tmpFR;
                
                % tuning curve for location
                loc_corrFR      = nan(numel(place.zCenters), numel(place.xCenters)); % z-axis: rows; x-axis: columns
                tmpFR           = empANOVA.m{strcmp(empANOVA.groupNames, 'P')}(:, 1);
                idx             = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'P')});
                loc_corrFR(idx) = tmpFR;
                loc_corrFR      = flipud(loc_corrFR); % flip up-down because index 1 shall be in the lower-left corner
                                
                %% surrogates
                
                % create surrogates or load previously saved surrogates
                resFile     = strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_surrogates.mat');
                if exist(resFile, 'file')
                    
                    % load previously saved results
                    tmp                     = load(resFile);
                    dir_Fsurro              = tmp.dir_Fsurro; % surrogate F values
                    all_dir_corrFRsurro     = tmp.all_dir_corrFRsurro; % surrogate tuning curves
                    loc_Fsurro              = tmp.loc_Fsurro;
                    all_loc_corrFRsurro     = tmp.all_loc_corrFRsurro;
                else
                    
                    % preallocate
                    dir_Fsurro              = nan(param.numSurrogates, 1);
                    all_dir_corrFRsurro     = nan(size(dir_corrFR, 1), param.numSurrogates);
                    loc_Fsurro              = nan(param.numSurrogates, 1);
                    all_loc_corrFRsurro     = nan(size(loc_corrFR, 1), size(loc_corrFR, 2), param.numSurrogates);
                    
                    % random shifts of firing rates
                    randShifts  = transpose(datasample(1:numel(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
                    
                    % loop through surrogates
                    parfor iSurro = 1:param.numSurrogates
                        
                        % create surrogate firing rate by circularly
                        % shifting firing rates relative to the behavioral
                        % data
                        FRsurro         = circshift(FR(bMask4Analysis), randShifts(iSurro, 1));
                        
                        %% surrogate ANOVA: direction and place on firing rate
                        
                        % compute ANOVAN
                        dt              = [];
                        dt.FR           = FRsurro; % you already masked the FR by the analysis period (see above)
                        dt.group        = [outBeh.yawBins(bMask4Analysis), outBeh.xzBins(bMask4Analysis)];
                        dt.groupNames   = {'D'; 'P'}; % predictors
                        dt.ANOVA        = ANOVA;
                        surroANOVA      = LK_TH_CompANOVA_210520(dt);
                        
                        % F values
                        dir_Fsurro(iSurro, 1)    = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'D'), strcmp(surroANOVA.tbl(1, :), 'F')));
                        loc_Fsurro(iSurro, 1)    = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'P'), strcmp(surroANOVA.tbl(1, :), 'F')));
                        
                        %% tuning curves
                        
                        % tuning curve for direction
                        dir_corrFRsurro                     = nan(numel(direc.angleCenters), 1);
                        tmpFR                               = surroANOVA.m{strcmp(surroANOVA.groupNames, 'D')}(:, 1);
                        idx                                 = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'D')});
                        dir_corrFRsurro(idx)                = tmpFR;
                        all_dir_corrFRsurro(:, iSurro)      = dir_corrFRsurro;
                        
                        % tuning curve for location
                        loc_corrFRsurro                     = nan(numel(place.zCenters), numel(place.xCenters));
                        tmpFR                               = surroANOVA.m{strcmp(surroANOVA.groupNames, 'P')}(:, 1);
                        idx                                 = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'P')});
                        loc_corrFRsurro(idx)                = tmpFR;
                        all_loc_corrFRsurro(:, :, iSurro)   = flipud(loc_corrFRsurro); % flip up-down because index 1 shall be in the lower-left corner                        
                    end
                    
                    % save surrogates
                    save(resFile, ...
                        'dir_Fsurro', 'all_dir_corrFRsurro', 'loc_Fsurro', 'all_loc_corrFRsurro');
                end
                
                %% statistical evaluation
                
                % rank of directional F-value
                dir_Frank 	= sum(dir_F > dir_Fsurro) / sum(~isnan(dir_Fsurro));
                
                % rank of locational F-value
                loc_Frank 	= sum(loc_F > loc_Fsurro) / sum(~isnan(loc_Fsurro));
                
                %% receptive bins for place (cf. Ekstrom et al., 2003)
                
                % location bins
                locBins     = (sum(loc_corrFR > all_loc_corrFRsurro, 3) ./ sum(~isnan(all_loc_corrFRsurro), 3)) > place.cutoffSigBin;
                
                %% collect information for this unit
                
                % original session index
                sessIdx     = split(sessions(iSess).name, '_');
                sessIdx     = str2double(sessIdx{2});
                
                % basics
                unitRes                 = [];
                unitRes.idx            	= [iSub, sessIdx, iWire, iClus];
                unitRes.bMask4Analysis  = bMask4Analysis;
                unitRes.wireRegion     	= wireRegion;
                unitRes.wireMNI        	= wireMNI;
                
                % direction
                unitRes.dir_F         	= dir_F;
                unitRes.dir_Fsurro    	= dir_Fsurro;
                unitRes.dir_Frank     	= dir_Frank;
                unitRes.dir_corrFR    	= dir_corrFR;
                
                % place
                unitRes.loc_F        	= loc_F;
                unitRes.loc_Fsurro    	= loc_Fsurro;
                unitRes.loc_Frank     	= loc_Frank;
                unitRes.loc_corrFR     	= loc_corrFR;
                unitRes.locBins        	= locBins;
                
                % save unit result
                save(strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_unitRes'), 'unitRes');
                
                % collapse across units
                allRes  = cat(1, allRes, unitRes);
                
                %% figure: direction
                
                % data for figure
                dt                  = [];
                dt.visible          = 'off';
                dt.figEx            = [2, 22, 8.5, 8]; % figure extension
                dt.circleEx         = [2.2, 0.4, 5.2, 5.2];
                dt.waveEx           = [1, 5.5, 2, 2];
                dt.direc            = direc;
                dt.FR               = dir_corrFR; % tuning curve
                dt.angLabels        = {'0°', '90°', '±180°', '-90°'};
                dt.bReverseX        = true; % whether to flip the x-axis
                dt.bReverseZ        = true; % whether to flip the z-axis
                dt.thisSpike        = thisSpike; % spike waveforms
                dt.t                = t; % sampling rate
                dt.numSpikes        = numSpikes;
                dt.bMask4Analysis   = bMask4Analysis;
                dt.figTitle         = {'Direction', ...
                    ['\itP\rm = ', num2str(1 - dir_Frank, '%.2f'), ' (', wireRegion, ')']};
                if (1 - dir_Frank) < 0.01
                    dt.figTitle     = {'Direction', ...
                        ['\itP\rm < 0.01 (', wireRegion, ')']};
                elseif (1 - dir_Frank) < 0.05 && round(1 - dir_Frank, 2) == 0.05
                    dt.figTitle     = {'Direction', ...
                        ['\itP\rm < 0.05 (', wireRegion, ')']};
                end
                
                % plot figure
                f = LK_TH_PlotFRDir_170520(dt);
                
                % save figure
                set(f, 'PaperPositionMode', 'auto');
                print(f, strcat(paths.save, subjects{iSub}, '_', sessions(iSess).name, '_', wires(iWire).name, '_clus', num2str(iClus), '_Dir'), '-dtiff', '-r150');
                close(f);
                
                %% figure: place
                
                % augmentation factor
                augFac          = 20;
                
                % resize x- and y-axes
                res_xCenters    = fliplr(movmean(linspace(min(place.xEdges), max(place.xEdges), augFac * numel(place.xCenters) + 1), 2, 'endpoints', 'discard'));
                res_zCenters    = fliplr(movmean(linspace(min(place.zEdges), max(place.zEdges), augFac * numel(place.zCenters) + 1), 2, 'endpoints', 'discard'));
                
                % resize the place map matrix and smooth
                resFR           = myresizem(loc_corrFR, augFac);
                dt              = [];
                dt.FR           = resFR;
                dt.smoothFac    = 21;
                outSmooth       = LK_TH_locSmooth_170520(dt);
                
                % resize location bins
                resLocBins      = myresizem(locBins, augFac);
                
                % data for figure
                dt                  = [];
                dt.visible          = 'off';
                dt.figEx            = [2, 22, 8.5, 8]; % figure extension
                dt.circleEx         = [2.6, 0.4, 5.5, 5.5];
                dt.waveEx           = [1, 5.5, 2, 2];
                dt.xCenters         = res_xCenters;
                dt.zCenters         = res_zCenters;
                dt.FR               = outSmooth.smFR; % firing rate
                dt.locBins          = resLocBins;
                dt.outBeh           = outBeh;
                dt.thisSpike        = thisSpike; % spike waveforms
                dt.t                = t; % sampling rate
                dt.numSpikes        = numSpikes;
                dt.bMask4Analysis   = bMask4Analysis;
                dt.figTitle         = {'Place', ...
                    ['\itP\rm = ', num2str(1 - loc_Frank, '%.2f'), ' (', wireRegion, ')']};
                if (1 - loc_Frank) < 0.01
                    dt.figTitle   	= {'Place', ...
                        ['\itP\rm < 0.01 (', wireRegion, ')']};
                elseif (1 - loc_Frank) < 0.05 && round(1 - loc_Frank, 2) == 0.05
                    dt.figTitle   	= {'Place', ...
                        ['\itP\rm < 0.05 (', wireRegion, ')']};
                end
                
                % plot figure
                f = LK_TH_plotFRLoc_170520(dt);
                
                % save figure
                set(f, 'PaperPositionMode', 'auto');
                print(f, strcat(paths.save, subjects{iSub}, '_', sessions(iSess).name, '_', wires(iWire).name, '_clus', num2str(iClus), '_Loc'), '-dtiff', '-r450');
                close(f);
                
                %% close all open figures
                close all;
                toc
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'), ...
    'allRes', 'allBeh', 'paths', ...
    'param', 'direc', 'place', 'subjects');

% all cells' indices
allSubjIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');
fprintf('Total number of cells: %d.\n', size(allSubjIdx, 1));

%% direction and place cells

% identify direction and place cells
bDirCell    = cell2mat({r.allRes.dir_Frank}') > 0.95;
bPlaceCell  = cell2mat({r.allRes.loc_Frank}') > 0.95;

% report
fprintf('\nNumber of direction cells: %d (P = %.3f).\n', sum(bDirCell), ...
    myBinomTest(sum(bDirCell), numel(bDirCell), 0.05));
fprintf('Number of place cells: %d (P = %.3f).\n', sum(bPlaceCell), ...
    myBinomTest(sum(bPlaceCell), numel(bPlaceCell), 0.05));
