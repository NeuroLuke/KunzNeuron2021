%==========================================================================
% This script examines whether cells change their firing rates during
% location-cued object recall.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath('E:\TreasureHunt\Functions\');
addpath('E:\OpenField\Functions\');

% variables
param                   = [];
param.clusterName       = 'Cluster4Analysis_LK_20200416_113912';
param.timeRes           = 0.25; % temporal resolution
param.trialTimeBorders 	= -1:param.timeRes:5; % relative to onset of location cue
param.trialTime         = movmean(param.trialTimeBorders, 2, 'endpoints', 'discard');
param.baselineBorders   = [-1, 0];
param.smoothKernel      = 'gaussian'; % smoothing kernel
param.smoothFac         = 1; % amount of smoothing (1 = no smoothing)
param.numSurrogates     = 1001; % number of surrogates
param.myRNG             = 4444;

% smoothing after averaging across trials
param.afterMeanSmoothKernel = 'gaussian';
param.afterMeanSmoothFac    = 1; % amount of smoothing (1 = no smoothing)

% paths
paths           = [];
paths.spike     = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.beh       = 'E:\TreasureHunt\Beh_210420\';
paths.save      = strcat('E:\TreasureHunt\ObjRecall_20200622\20210418\');
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% get subjects
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

%% load previous results

% load previous ELRPD cell results
ELRPD   = load('E:\TreasureHunt\ELRPDCellAnalysis_220420\20200722_Dir30_Loc6x6_EgoLocDir30\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', 'allRes');

%% preallocations

% main results for each cell
allRes   	= [];

% for sanity checking
exByWC      = []; % excluded by wave-clus
exByELRPD   = []; % excluded due to ELRPD-cell analysis

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % available sessions
    sessions 	= dir(strcat(paths.spike, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        %% behavioral data
        
        % original session index
        sessIdx     = split(sessions(iSess).name, '_');
        sessIdx     = str2double(sessIdx{2});
        
        % report
        fprintf('\nSubject: %s, session: %s ...\n', subjects{iSub}, sessions(iSess).name);
        
        % load behavioral data
        tI          = load(strcat(paths.beh, subjects{iSub}, '\', sessions(iSess).name, '\trialInfo.mat'));
        trialInfo   = tI.trialInfo;
        
        % object-recall performance
        objRecallPerf   = trialInfo.CORTANA_RESPONSE;
        
        %% wires to investigate
        
        % get available microwires
        wires       = dir(fullfile(sessions(iSess).folder, sessions(iSess).name, 'chan*'));
        tmp         = split(transpose({wires.name}), 'n');
        [~, I]      = sort(cellfun(@str2num, tmp(:, 2)));
        wires       = wires(I);
        
        %% loop through wires
        for iWire = 1:size(wires, 1)
            
            % report
            fprintf('\tWire: %s.\n', wires(iWire).name);
            
            %% get spike times in behavioral time
            
            % load wave-clus output (for spike depiction etc.)
            wirePath    = strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep);
            try
                t = load(strcat(wirePath, 'times_datacut.mat'));
            catch
                exByWC  = cat(1, exByWC, [iSub, sessIdx, iWire]); % bookkeeping
                fprintf('No wave-clus for this wire.\n');
                continue;
            end
            
            % load decision whether to use clusters
            c4a     = load(strcat(wirePath, param.clusterName, '.mat'));
            
            % load behavioral times of spikes
            ccb     = load(strcat(wirePath, 'cluster_class_behtime.mat'));
            
            %% loop through clusters
            for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
                
                % report
                fprintf('\t\tCluster: %d.\n', iClus);
                
                % continue if this unit was not included in the ELRPD cell
                % analysis
                if ~any(all([iSub, sessIdx, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2))
                    exByELRPD   = cat(1, exByELRPD, [iSub, sessIdx, iWire, iClus]);
                    fprintf('This cluster was not included in the ELRPD-cell analysis.\n');
                    continue;
                end
                
                % data for this cluster
                thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster index, microwire time (msec), behavioral time (sec)
                
                %% firing rates during the baseline period
                
                % onsets, ends, and durations
                onsEnds     = trialInfo.RECORDING_STARTED + param.baselineBorders;
                durations   = diff(onsEnds, [], 2);
                
                % number of spikes and firing rates
                numSpikes   = nan(size(onsEnds, 1), 1);
                for iOns = 1:size(onsEnds, 1)
                    numSpikes(iOns, :)  = histcounts(thisCluster(:, 3), onsEnds(iOns, :));
                end
                preRecFR    = numSpikes ./ durations;
                
                %% firing rates during object recall
                
                % time bins
                trialOns            = trialInfo.RECORDING_STARTED;
                trialTimeBorders    = trialOns + param.trialTimeBorders; % time within each trial
                
                % number of spikes and durations
                numSpikes   = nan(size(trialTimeBorders, 1), size(trialTimeBorders, 2) - 1);
                for iTrial = 1:size(trialTimeBorders, 1)
                    numSpikes(iTrial, :)    = histcounts(thisCluster(:, 3), trialTimeBorders(iTrial, :));
                end
                durations       = diff(trialTimeBorders, [], 2);
                
                % firing rate
                objRecTrialFR   = numSpikes ./ durations;
                
                %% collect information for this unit
                
                % basics
                unitRes                 = [];
                unitRes.idx             = [iSub, sessIdx, iWire, iClus];
                
                % performance
                unitRes.objRecallPerf   = objRecallPerf;
                
                % firing rates during baseline
                unitRes.preRecFR        = preRecFR;
                
                % firing rates during recall
                unitRes.objRecTrialFR   = objRecTrialFR;
                
                % collapse across units
                allRes  = cat(1, allRes, unitRes);                
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));
fprintf('Total number of cells: %d.\n', size(r.allRes, 1));

% ELRPD cells
bELRPDCell  = cell2mat({r.ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells (from ELRPD-cell analysis): %d.\n', sum(bELRPDCell));

%% direction cells and place cells

% load results from previous place x direction analysis
PxD             = load('E:\TreasureHunt\PlaceDirAnalysis_220520\20200623_Dir30_Loc6x6\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', 'allRes');
bDirCell     	= cell2mat({PxD.allRes.dir_Frank}') > 0.95; % direction cells
bPlaceCell    	= cell2mat({PxD.allRes.loc_Frank}') > 0.95; % place-like cells
fprintf('Number of direction cells (from place*direction analysis): %d.\n', sum(bDirCell));
fprintf('Number of place cells (from place*direction analysis): %d.\n', sum(bPlaceCell));

% spatial cells (i.e., either ELRPD or direction or place cells)
bSpatialCell   	= bELRPDCell | bDirCell | bPlaceCell;
fprintf('Number of spatial cells: %d.\n', sum(bSpatialCell));

%% investigate time-resolved FRs

% average firing-rate timecourse for each unit, separately for good vs. bad
% trials
allMeanFR_good  = nan(size(r.allRes, 1), numel(r.param.trialTime)); % good trials
allMeanFR_bad   = nan(size(r.allRes, 1), numel(r.param.trialTime)); % bad trials
for iCell = 1:size(r.allRes, 1)
    
    % memory performance
    bGoodMemory = strcmp(r.allRes(iCell).objRecallPerf, 'CORRECT');
    
    % smoothing
    smFR    = smoothdata(r.allRes(iCell).objRecTrialFR, 2, r.param.smoothKernel, r.param.smoothFac); % smooth along second dimension
    
    % baseline correction
    logIdx  = r.param.trialTime >= r.param.baselineBorders(1) & r.param.trialTime <= r.param.baselineBorders(2);
    smFRbc 	= smFR - mean(smFR(:, logIdx), 2);
    
    % collect across cells
    allMeanFR_good(iCell, :)    = mean(smFRbc(bGoodMemory, :), 1);
    allMeanFR_bad(iCell, :)     = mean(smFRbc(~bGoodMemory, :), 1);
end

% smoothing after averaging
allMeanFR_good  = smoothdata(allMeanFR_good, 2, r.param.afterMeanSmoothKernel, r.param.afterMeanSmoothFac); % smooth along second dimension
allMeanFR_bad   = smoothdata(allMeanFR_bad, 2, r.param.afterMeanSmoothKernel, r.param.afterMeanSmoothFac);

%% fieldtrip

% add fieldtrip to path
addpath('E:\fieldtrip\fieldtrip-20201113\');
ft_defaults;
ft_warning('off');

%% object recall, good vs. bad, cluster-based permutation test and figure

% report
fprintf('\n==== Testing firing rates during good vs. bad object recall, separately for different cell groups.\n');

% data groups and settings
groups  = {'spatialCells'; 'nonSpatialCells'; 'ELRPDCells'};
myYLim  = [-0.3, 0.9]; % y limits for figure

% loop through groups
for iGroup = 1:numel(groups)
    
    % reset rng
    rng(r.param.myRNG);
    
    % select specific cells
    fprintf('\n... %s:\n', groups{iGroup});
    if strcmp(groups{iGroup}, 'spatialCells')
        bCellMask   = bSpatialCell;
    elseif strcmp(groups{iGroup}, 'nonSpatialCells')
        bCellMask   = ~bSpatialCell;
    elseif strcmp(groups{iGroup}, 'ELRPDCells')
        bCellMask   = bELRPDCell;
    end
    
    % select data from specific cell group
    dataGood    = allMeanFR_good(bCellMask, :);
    dataBad     = allMeanFR_bad(bCellMask, :);
        
    %% cluster-based permutation testing using fieldtrip, good > bad
    
    % re-organize data for fieldtrip
    timelock1   = cell(1, size(dataGood, 1));
    timelock2  	= cell(1, size(dataBad, 1));
    for iCell = 1:size(dataGood, 1)
        
        % good trials
        timelock1{iCell}.avg       = dataGood(iCell, :);
        timelock1{iCell}.label     = {'neuron'};
        timelock1{iCell}.time      = r.param.trialTime;
        
        % bad trials
        timelock2{iCell}.avg       = dataBad(iCell, :);
        timelock2{iCell}.label     = {'neuron'};
        timelock2{iCell}.time      = r.param.trialTime;
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = 'all';
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'depsamplesT';
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05; % first-level alpha
    cfg.clusterstatistic    = 'maxsum';
    cfg.neighbours          = [];
    cfg.tail                = 1; % H1: stronger activation during good trials
    cfg.alpha               = 0.05; % second-level alpha
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = r.param.numSurrogates; % number of surrogates
    numCells                = size(dataGood, 1);
    design                  = zeros(2, numCells * 2);
    design(1, :)            = [1:numCells, 1:numCells]; % cell indices
    design(2, :)            = [ones(1, numCells), ones(1, numCells) * 2]; % group indices
    cfg.design              = design; % design matrix
    cfg.uvar                = 1;
    cfg.ivar                = 2;
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    
    % report p-value of largest positive cluster
    if isfield(outFT, 'posclusters') && ~isempty(outFT.posclusters)
        fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT.posclusters.prob})));
    end
        
    %% cluster-based permutation testing, good > 0
    
    % test for significant activation above zero
    c                   = [];
    c.mat              	= dataGood;
    c.alpha          	= 0.05;
    c.direction     	= 'pos'; % H1: positive activation during good trials
    c.numSurrogates    	= r.param.numSurrogates; % number of surrogates
    outCBPT             = LK_PermutationTest_OneSampleAgainstZero_20200708(c);
        
    %% figure
    
    % mean and standard error
    m   = {mean(dataBad, 1), mean(dataGood, 1)}; % bad, good
    sem = {std(dataBad, 1) ./ sqrt(sum(~isnan(dataBad), 1)), ...
        std(dataGood, 1) ./ sqrt(sum(~isnan(dataGood), 1))}; % bad, good
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{1} - sem{1}, fliplr(m{1} + sem{1})], [1, 0, 0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{2} - sem{2}, fliplr(m{2} + sem{2})], [0, 0.5, 0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(r.param.trialTime, m{1}, '-', ...
        'Color', [1, 0, 0], 'LineWidth', 2);
    plot(r.param.trialTime, m{2}, '-', ...
        'Color', [0, 0.5, 0], 'LineWidth', 2);
    set(gca, ...
        'YLim', myYLim);
    % significance info
    LK_SigLine(r.param.trialTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT.mask); % good > bad
    LK_SigLine(r.param.trialTime, [max(myYLim) - range(myYLim) * 0.04; max(myYLim) - range(myYLim) * 0.065], outCBPT.logIdxSigClus, [0.5, 0.5, 0.5]); % good > 0
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(param.trialTimeBorders), max(param.trialTimeBorders)], 'xtick', min(param.trialTimeBorders):max(param.trialTimeBorders), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(paths.save, groups{iGroup}, '_objRecall_trFR_goodVSbad_20201115'), '-dtiff', '-r300');
end

%% cluster-based permutation test for interaction (cell type x performance) using fieldtrip

% report
fprintf('\n==== Testing for an interaction between cell type and performance, using cluster-based permutation testing across the entire time window.\n');

% groups of cells that you want to contrast
groups  = {'spatialCells', 'nonSpatialCells'; ... % group A, group B
    'ELRPDCells', 'nonSpatialCells'};

% loop through groups
for iGroup = 1:size(groups, 1)
    
    % report
    fprintf('\n... Interaction good vs. bad trials and %s vs %s.\n\n', groups{iGroup, 1}, groups{iGroup, 2});
    
    % difference between conditions, separately for the two groups
    if strcmp(groups{iGroup, 1}, 'ELRPDCells')
        groupA_goodVSbad 	= allMeanFR_good(bELRPDCell, :) - allMeanFR_bad(bELRPDCell, :); % ELRPD cells, good - bad
    elseif strcmp(groups{iGroup, 1}, 'spatialCells')
        groupA_goodVSbad    = allMeanFR_good(bSpatialCell, :) - allMeanFR_bad(bSpatialCell, :); % spatial cells, good - bad
    end
    if strcmp(groups{iGroup, 2}, 'nonSpatialCells')
        groupB_goodVSbad    = allMeanFR_good(~bSpatialCell, :) - allMeanFR_bad(~bSpatialCell, :); % non-spatial cells, good - bad
    end
    
    % re-organize data for fieldtrip
    timelock1       = [];
    timelock1.label = {'neuron'};
    for iCell = 1:size(groupA_goodVSbad, 1)
        timelock1.time{1, iCell}    = r.param.trialTime;
        timelock1.trial{1, iCell}   = groupA_goodVSbad(iCell, :);
    end
    timelock2       = [];
    timelock2.label = {'neuron'};
    for iCell = 1:size(groupB_goodVSbad, 1)
        timelock2.time{1, iCell}    = r.param.trialTime;
        timelock2.trial{1, iCell}   = groupB_goodVSbad(iCell, :);
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = 'all';
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'indepsamplesT'; % first-level test
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05; % first-level alpha
    cfg.clusterstatistic    = 'maxsum'; % test statistic
    cfg.neighbours          = [];
    cfg.tail                = 1; % H1: stronger activation during good trials for ELRPD cells and spatial cells
    cfg.alpha               = 0.05; % second-level alpha
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = r.param.numSurrogates;
    cfg.design              = [ones(size(timelock1.trial)), 2 .* ones(size(timelock2.trial))]; % design matrix
    cfg.ivar                = 1;
    
    % fieldtrip estimation
    outFT_inter = ft_timelockstatistics(cfg, timelock1, timelock2);
    
    % report p-value of largest positive cluster
    if isfield(outFT_inter, 'posclusters') && ~isempty(outFT_inter.posclusters)
        fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT_inter.posclusters.prob})));
        fprintf('The significant window lasts from %.3f to %.3f.\n', min(r.param.trialTime(outFT_inter.mask)), max(r.param.trialTime(outFT_inter.mask)));
    end
    
    %% figure for interaction effect
    
    % mean and standard error
    m   = {mean(groupB_goodVSbad, 1), mean(groupA_goodVSbad, 1)}; % difference #2, difference #1
    sem = {std(groupB_goodVSbad, 1) ./ sqrt(sum(~isnan(groupB_goodVSbad), 1)), ...
        std(groupA_goodVSbad, 1) ./ sqrt(sum(~isnan(groupA_goodVSbad), 1))}; % difference #2, difference #1
    
    % colors
    myColor1    = [0.7, 0.7, 0.7];
    myColor2    = [0.4, 0.4, 0.4];
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{1} - sem{1}, fliplr(m{1} + sem{1})], myColor1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{2} - sem{2}, fliplr(m{2} + sem{2})], myColor2, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(r.param.trialTime, m{1}, '-', ...
        'Color', myColor1, 'LineWidth', 2);
    plot(r.param.trialTime, m{2}, '-', ...
        'Color', myColor2, 'LineWidth', 2);
    set(gca, ...
        'YLim', myYLim);
    % significance info
    LK_SigLine(r.param.trialTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT_inter.mask);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(param.trialTimeBorders), max(param.trialTimeBorders)], 'xtick', min(param.trialTimeBorders):max(param.trialTimeBorders), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(paths.save, 'ObjRecall_InteractionEffectCBPT_goodVSbad_', groups{iGroup, 1}, 'VS', groups{iGroup, 2}), '-dtiff', '-r600');
end
