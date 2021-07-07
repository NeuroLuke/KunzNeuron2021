%==========================================================================
% This script analyzes egocentric bearing cells (= ELRPD cells).
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
param.newTimeRes        = 0.1; % temporal resolution
param.naviCutoff        = 0.1; % speed cutoff
param.naviSmooth        = 2;
param.rotationCutoff    = 0.01; % angular speed cutoff
param.rotationSmooth    = 2;
param.arenaXLim         = [420, 320]; % west, east
param.arenaZLim         = [410, 310]; % south, north
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
place.locRes            = 6; % spatial resolution (cf. Tsitsiklis et al., 2020)
place.idxTemplate       = flipud(reshape(1:place.locRes ^ 2, place.locRes, place.locRes)); % template for location indexing
place.xEdges            = linspace(430, 310, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.zEdges            = linspace(420, 300, place.locRes + 1);
place.zCenters          = movmean(place.zEdges, 2, 'endpoints', 'discard');
place.cutoffSigBin      = 0.95;
%- for egocentric local directions
locDir                  = [];
locDir.locRes           = 12; % spatial resolution
locDir.xEdges           = linspace(430, 310, locDir.locRes + 1);
locDir.xCenters         = movmean(locDir.xEdges, 2, 'endpoints', 'discard');
locDir.zEdges           = linspace(420, 300, locDir.locRes + 1);
locDir.zCenters         = movmean(locDir.zEdges, 2, 'endpoints', 'discard');
locDir.angularRes       = 30; % angular resolution
locDir.angleEdges       = deg2rad(-180:locDir.angularRes:180);
locDir.angleCenters     = movmean(locDir.angleEdges, 2, 'endpoints', 'discard');
%- local reference points
tmpX                    = meshgrid(locDir.xCenters);
tmpZ                    = transpose(meshgrid(locDir.zCenters));
locDir.LRPs             = [tmpX(:), tmpZ(:)];
% remove LRPs that are too far from center
D_LRPs2Ctr              = pdist2(locDir.LRPs, param.arenaCtr);
locDir.LRPs2Exclude     = D_LRPs2Ctr > 60;
%- for location-specific directional tuning
locDir.areaAroundLRP    = param.arenaRadius * 2/3; % 1/3 of the diameter
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sums of squares

% paths
paths           = [];
paths.info      = 'E:\TreasureHunt\SubjectInformation_170818\';
paths.spike     = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.beh       = 'E:\TreasureHunt\Beh_210420\';
paths.save      = strcat('E:\TreasureHunt\ELRPDCellAnalysis_220420\20200722_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.zCenters)), ...
    '_EgoLocDir', num2str(locDir.angularRes), '\', ...
    'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotationIfAbove', num2str(param.rotationCutoff), '_Smooth', num2str(param.naviSmooth), '\');
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
    
    % available sessions
    sessions    = dir(fullfile(paths.spike, subjects{iSub}, 'session_*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % report
        fprintf('\n\nSUBJECT: %s, SESSION: %s.\n', subjects{iSub}, sessions(iSess).name);
        rng(param.myRNG);
        
        % original session index
        sessIdx = split(sessions(iSess).name, '_');
        sessIdx = str2double(sessIdx{2});
        
        % subject-specific save path
        subjSavePath    = fullfile(paths.save, subjects{iSub}, sessions(iSess).name, '\');
        if ~exist(subjSavePath, 'dir')
            mkdir(subjSavePath); % create subject-specific folder
        end
        
        %% behavioral data
        
        % process behavioral data
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
            outBeh.xBins    = numel(place.xEdges) - discretize(outBeh.xyz(:, 1), fliplr(place.xEdges)); % high x-values get low bin indices
            outBeh.zBins    = numel(place.zEdges) - discretize(outBeh.xyz(:, 3), fliplr(place.zEdges)); % high z-values get low bin indices
            outBeh.xzBins   = (outBeh.xBins - 1) .* numel(place.zCenters) + outBeh.zBins;
            
            % egocentric local directions for all LRPs
            egoLocDir       = cell(size(locDir.LRPs, 1), 1);
            egoLocDirBins   = cell(size(locDir.LRPs, 1), 1);
            parfor iLRP = 1:size(locDir.LRPs, 1)
                                
                % calculate egocentric directions
                alloLocDir          = atan2(locDir.LRPs(iLRP, 2) - outBeh.xyz(:, 3), locDir.LRPs(iLRP, 1) - outBeh.xyz(:, 1)); % atan2(y, x)
                egoLocDir{iLRP, 1}  = angdiff(outBeh.yaws, alloLocDir); % angdiff(x, y) returns the angular difference delta = y - x
                
                % discretize egocentric directions
                egoLocDirBins{iLRP, 1}  = discretize(egoLocDir{iLRP, 1}, locDir.angleEdges);
            end
            
            % add to other behavioral information
            outBeh.egoLocDir        = egoLocDir;
            outBeh.egoLocDirBins    = egoLocDirBins;
            
            % save output
            save(behFile, 'outBeh');
        end
        
        % collect behavioral output for all subjects and sessions
        allBeh  = cat(1, allBeh, outBeh);
        
        %% mask for analysis timepoints
        
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
        
        % boolean that encodes whether specific time bins shall be included
        % in the analysis
        bMask4Analysis  = strcmp(outBeh.behInfo(:, 10), 'navigation') & (outBeh.bNavi | outBeh.bRotation);
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
                exByWC  = cat(1, exByWC, [iSub, sessIdx, iWire]); % bookkeeping
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
                    exByVisInsp = cat(1, exByVisInsp, [iSub, sessIdx, iWire, iClus]); % bookkeeping
                    continue;
                end
                
                % data from this cluster
                thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
                thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
                
                %% number of spikes and FR per time bin
                
                % number of spikes per timepoint
                numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time-bin
                
                % firing rate per timepoint
                FR  = numSpikes ./ outBeh.durations;
                
                % if the general firing rate is too low (<0.1 Hz), continue
                if (sum(numSpikes(bMask4Analysis)) / sum(outBeh.durations(bMask4Analysis))) < 0.1
                    fprintf('\t\t- The mean firing rate of this cluster is less than 0.1 Hz, thus skipping.\n');
                    exByNumSpikes   = cat(1, exByNumSpikes, [iSub, sessIdx, iWire, iClus]); % bookkeeping
                    continue;
                end
                                
                %% egocentric bearing tuning towards each LRP
                
                % for each local reference point, estimate modulation of
                % firing rate by egocentric direction towards this LRP
                resFile     = strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_allLocDir_emp.mat');
                if exist(resFile, 'file')
                    
                    % load previously saved results
                    tmp                 = load(resFile);
                    all_locdir_F        = tmp.all_locdir_F; % F values
                    all_locdir_corrFR   = tmp.all_locdir_corrFR; % tuning curves
                else
                    
                    % preallocate
                    all_locdir_F        = nan(size(locDir.LRPs, 1), 1);
                    all_locdir_corrFR   = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters));
                    for iLRP = 1:size(locDir.LRPs, 1)
                        
                        % compute ANOVAN
                        dt              = [];
                        dt.FR           = FR(bMask4Analysis);
                        dt.group        = [outBeh.yawBins(bMask4Analysis), outBeh.xzBins(bMask4Analysis), outBeh.egoLocDirBins{iLRP, 1}(bMask4Analysis)];
                        dt.groupNames   = {'D'; 'P'; 'ELRPD'}; % predictors
                        dt.ANOVA        = ANOVA;
                        empANOVA        = LK_TH_CompANOVA_210520(dt);
                        
                        % F value
                        all_locdir_F(iLRP, 1)   = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'ELRPD'), strcmp(empANOVA.tbl(1, :), 'F')));
                        
                        % tuning curve for ELRPD
                        tmpFR                           = empANOVA.m{strcmp(empANOVA.groupNames, 'ELRPD')}(:, 1);
                        idx                             = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'ELRPD')});
                        all_locdir_corrFR(iLRP, idx)    = tmpFR;
                    end
                    
                    % save result
                    save(resFile, ...
                        'all_locdir_F', 'all_locdir_corrFR');
                end
                
                %% surrogates for egocentric bearing tuning towards each LRP
                
                % random numbers for circular shift of firing rates
                randShifts  = transpose(datasample(1:numel(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
                
                % create surrogates or load previously saved surrogates
                resFile     = strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_allLocDir_surrogates.mat');
                if exist(resFile, 'file')
                    
                    % load previously saved results
                    tmp                     = load(resFile);
                    all_locdir_Fsurro       = tmp.all_locdir_Fsurro;
                    all_locdir_corrFRsurro  = tmp.all_locdir_corrFRsurro;
                else
                    
                    % preallocate
                    all_locdir_Fsurro      	= nan(size(locDir.LRPs, 1), param.numSurrogates); % e.g., size = 144 x 101
                    all_locdir_corrFRsurro  = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters), param.numSurrogates); % e.g., size = 144 x 12 x 101
                     
                    % loop through surrogates
                    parfor iSurro = 1:param.numSurrogates
                        
                        % create surrogate firing rate by circularly
                        % shifting firing rates relative to the behavioral
                        % data
                        FRsurro = circshift(FR(bMask4Analysis), randShifts(iSurro, 1));
                                                
                        %% surrogate ANOVA for modulation of firing rate by egocentric direction towards local reference points
                        
                        % for each local reference point, estimate
                        % modulation of firing rate by egocentric direction
                        % towards this LRP
                        tmp_locdir_Fsurro       = nan(size(locDir.LRPs, 1), 1);
                        tmp_locdir_corrFRsurro  = nan(size(locDir.LRPs, 1), numel(locDir.angleCenters));
                        for iLRP = 1:size(locDir.LRPs, 1)
                            
                            % compute ANOVAN
                            dt              = [];
                            dt.FR           = FRsurro; % you already masked the FR by the analysis period (see above)
                            dt.group        = [outBeh.yawBins(bMask4Analysis), outBeh.xzBins(bMask4Analysis), outBeh.egoLocDirBins{iLRP, 1}(bMask4Analysis)];
                            dt.groupNames   = {'D'; 'P'; 'ELRPD'};
                            dt.ANOVA        = ANOVA;
                            surroANOVA      = LK_TH_CompANOVA_210520(dt);
                            
                            % F value for ELRPD
                            tmp_locdir_Fsurro(iLRP, 1)  = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'ELRPD'), strcmp(surroANOVA.tbl(1, :), 'F')));
                            
                            % tuning curve for ELRPD
                            tmpFR                               = surroANOVA.m{strcmp(surroANOVA.groupNames, 'ELRPD')}(:, 1);
                            idx                                 = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'ELRPD')});
                            tmp_locdir_corrFRsurro(iLRP, idx)   = tmpFR;
                        end
                        
                        % collect across surrogate rounds
                        all_locdir_Fsurro(:, iSurro)          	= tmp_locdir_Fsurro;
                        all_locdir_corrFRsurro(:, :, iSurro)    = tmp_locdir_corrFRsurro;
                    end
                    
                    % save surrogates
                    save(resFile, ...
                        'all_locdir_Fsurro', 'all_locdir_corrFRsurro');
                end
                
                %% statistical evaluation of ELRPD tuning
                
                % rank of ELRPD F value for each LRP
                all_locdir_Frank        = sum(all_locdir_F > all_locdir_Fsurro, 2) ./ sum(~isnan(all_locdir_Fsurro), 2);
                
                % test whether the modulation map is significant using
                % cluster-based permutation testing (CBPT)
                dt                      = [];
                dt.locDir               = locDir;
                dt.all_locdir_F         = all_locdir_F;
                dt.all_locdir_Fsurro    = all_locdir_Fsurro;
                dt.all_locdir_Frank     = all_locdir_Frank;
                dt.numSurrogates        = param.numSurrogates;
                locdir_CBPT             = LK_TH_locDir_CBPT_280520(dt);
                
                %% location-specific allocentric direction tuning at each LRP (+ surrounding area)
                
                % create or load output
                resFile = strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_dir_corrFR_perLRP.mat');
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
                        
                        % select time points when subject was in the
                        % vicinity of this LRP
                        bMask4AnalysisLRP   = bMask4Analysis & ...
                            outBeh.xyz(:, 1) >= (locDir.LRPs(iLRP, 1) - locDir.areaAroundLRP) & ...
                            outBeh.xyz(:, 1) <= (locDir.LRPs(iLRP, 1) + locDir.areaAroundLRP) & ...
                            outBeh.xyz(:, 3) >= (locDir.LRPs(iLRP, 2) - locDir.areaAroundLRP) & ...
                            outBeh.xyz(:, 3) <= (locDir.LRPs(iLRP, 2) + locDir.areaAroundLRP);
                        
                        % compute ANOVAN for this data part
                        try
                            dt              = [];
                            dt.FR           = FR(bMask4AnalysisLRP);
                            dt.group        = [outBeh.yawBins(bMask4AnalysisLRP), outBeh.xzBins(bMask4AnalysisLRP)];
                            dt.groupNames   = {'D'; 'P'}; % predictors
                            dt.ANOVA        = ANOVA;
                            empANOVA        = LK_TH_CompANOVA_210520(dt);
                            
                            % output F values
                            dir_F_perLRP(iLRP, 1)           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'D'), strcmp(empANOVA.tbl(1, :), 'F')));
                            
                            % tuning curve for direction
                            tmpFR                           = empANOVA.m{strcmp(empANOVA.groupNames, 'D')}(:, 1);
                            idx                             = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')});
                            dir_corrFR_perLRP(iLRP, idx)    = tmpFR;
                        catch
                            fprintf('\t\t--- Location-specific analysis of direction tuning was not possible for LRP (%d/%d).\n', ...
                                locDir.LRPs(iLRP, 1), locDir.LRPs(iLRP, 2));
                        end
                    end
                    
                    % save output
                    save(resFile, 'dir_F_perLRP', 'dir_corrFR_perLRP');
                end
                
                %% collect information for this unit
                                
                % basics
                unitRes                             = [];
                unitRes.idx                         = [iSub, sessIdx, iWire, iClus];
                unitRes.bMask4Analysis              = bMask4Analysis;
                unitRes.wireRegion                  = wireRegion;
                unitRes.wireMNI                     = wireMNI;
                
                % all ELRPD tunings
                unitRes.all_locdir_F                = all_locdir_F;
                unitRes.all_locdir_Fsurro           = all_locdir_Fsurro;
                unitRes.all_locdir_Frank            = all_locdir_Frank;
                unitRes.all_locdir_corrFR        	= all_locdir_corrFR;
                unitRes.locdir_clusRank             = locdir_CBPT.locdir_clusRank;
                unitRes.locdir_maxSumStat_emp       = locdir_CBPT.maxSumStat_emp;
                unitRes.locdir_largestCluster       = locdir_CBPT.largestCluster;
                unitRes.locdir_COMxz                = locdir_CBPT.COMxz;
                unitRes.locdir_COMxzClosestLRP      = locdir_CBPT.COMxzClosestLRP;
                
                % direction, separately for each LRP
                unitRes.dir_F_perLRP                = dir_F_perLRP;
                unitRes.dir_corrFR_perLRP           = dir_corrFR_perLRP;
                
                % save results of this unit
                save(strcat(subjSavePath, wires(iWire).name, '_clus', num2str(iClus), '_unitRes'), 'unitRes');
                
                % collapse across units
                allRes  = cat(1, allRes, unitRes);
                
                %% figure: ELRPD
                
                % data for figure
                dt                  	= [];
                dt.visible            	= 'off';
                dt.param             	= param;
                dt.locDir             	= locDir;
                dt.locdir_CBPT        	= locdir_CBPT;
                dt.all_locdir_corrFR   	= all_locdir_corrFR;
                dt.direc              	= direc;
                dt.dir_corrFR_perLRP  	= dir_corrFR_perLRP;
                dt.thisSpike           	= thisSpike; % spike waveforms
                dt.t.par.sr            	= t.par.sr; % sampling rate
                dt.nspk              	= sum(numSpikes(bMask4Analysis));
                dt.figTitle         	= {'Reference field and point', ...
                    ['\itP\rm = ', num2str(1 - locdir_CBPT.locdir_clusRank, '%.2f'), ' (', wireRegion, ')']};
                if (1 - locdir_CBPT.locdir_clusRank) < 0.05 && round(1 - locdir_CBPT.locdir_clusRank, 2) == 0.05
                    dt.figTitle        	= {'Reference field and point', ...
                        ['\itP\rm < 0.05 (', wireRegion, ')']};
                elseif (1 - locdir_CBPT.locdir_clusRank) < 0.01
                    dt.figTitle        	= {'Reference field and point', ...
                        ['\itP\rm < 0.01 (', wireRegion, ')']};
                end
                
                % plot figure
                f = LK_TH_PlotELRPDModulation_20210515(dt);
                
                % save figure
                set(f, 'PaperPositionMode', 'auto');
                print(f, strcat(paths.save, subjects{iSub}, '_', sessions(iSess).name, '_', wires(iWire).name, '_clus', num2str(iClus), '_ELRPD'), '-dtiff', '-r450');
                close(f);
                
                %% close all open figures
                close all;
                toc                
            end
        end
    end
end

%% save all results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'), ...
    'allRes', 'allBeh', 'paths', ...
    'param', 'direc', 'place', 'locDir', 'subjects');

% all cell indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');

%% load results from direction and place analysis

% load data
PxD             = load('E:\TreasureHunt\PlaceDirAnalysis_220520\20200623_Dir30_Loc6x6\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', ...
    'allRes');
PxD_bDirCell    = cell2mat({PxD.allRes.dir_Frank}') > 0.95; % direction cells
PxD_bPlaceCell  = cell2mat({PxD.allRes.loc_Frank}') > 0.95; % place-like cells

% report
fprintf('\nNumber of direction cells from the PxD analysis: %d (P = %.3f).\n', sum(PxD_bDirCell), ...
    myBinomTest(sum(PxD_bDirCell), numel(PxD_bDirCell), 0.05));
fprintf('Number of place-like cells from the PxD analysis: %d (P = %.3f).\n', sum(PxD_bPlaceCell), ...
    myBinomTest(sum(PxD_bPlaceCell), numel(PxD_bPlaceCell), 0.05));

%% ELRPD cells

% ELRPD cells
all_locdir_clusRank     = cell2mat({r.allRes.locdir_clusRank}');
bELRPDCell              = all_locdir_clusRank > 0.95; % >0.95 corresponds to P<0.05
fprintf('\nNumber of ELRPD cells: %d (P = %.3f).\n', sum(bELRPDCell), ...
    myBinomTest(sum(bELRPDCell), numel(bELRPDCell), 0.05));

% ELRPD cells that are not also direction cells or place-like cells
bPureELRPDCell 	= bELRPDCell & ~PxD_bDirCell & ~PxD_bPlaceCell;
fprintf('Number of pure ELRPD cells: %d (P = %.3f).\n', sum(bPureELRPDCell), ...
    myBinomTest(sum(bPureELRPDCell), numel(bPureELRPDCell), 0.05));

%% evaluate ELRPD cells

% characterize ELRPD cells
dt              = [];
dt.r            = r;
dt.bELRPDCell   = bELRPDCell;
dt.bDirCell     = PxD_bDirCell;
dt.bPlaceCell   = PxD_bPlaceCell;
LK_TH_EvalELRPDCells_20210515(dt);

%% region-wise distribution of cells

% report
fprintf('\nEvaluation of region-wise distribution of cells.\n');

% unique regions and number of cells per unique region
allRegions      = {r.allRes.wireRegion}';
uniqueRegions   = unique(allRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% loop through regions
for iReg = 1:size(uniqueRegions, 1)
    
    % identify cells from this region
    fprintf('\nRegion: %s.\n', uniqueRegions{iReg});
    bThisReg                = strcmp(allRegions, uniqueRegions{iReg});
    numCellsPerReg(iReg, 1) = sum(bThisReg);
    
    % ELRPD cells
    fprintf('ELRPD cells: %d (of %d, P = %.3f).\n', ...
        sum(bELRPDCell & bThisReg), sum(bThisReg), ...
        myBinomTest(sum(bELRPDCell & bThisReg), sum(bThisReg), 0.05)); 
end

% restrict unique regions to those with enough units
uniqueRegions   = uniqueRegions(numCellsPerReg >= 30);

%% bar plot for percentage of cells per region

% ELRPD cells
dt                  = [];
dt.figEx            = [5, 2, 8, 8];
dt.uniqueRegions    = uniqueRegions;
dt.allRegions       = allRegions;
dt.bCell            = bELRPDCell;
dt.ylabel           = 'Egocentric bearing cells (%)';
[f, percPerReg]     = LK_PlotCellsPerRegion(dt);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(r.paths.save, 'ELRPDCell_PerRegion_20210409'), '-dtiff', '-r300');

%% influence of recording site

% analyze effect of epileptic region and hemisphere
dt           	= [];
dt.r            = r;
dt.bELRPDCell   = bELRPDCell;
dt.allRegions   = allRegions;
LK_TH_EvalRecSites_20210515(dt);
