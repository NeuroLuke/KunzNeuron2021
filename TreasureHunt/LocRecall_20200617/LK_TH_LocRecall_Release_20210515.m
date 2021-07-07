%==========================================================================
% This script examines whether cells change their firing rates during 
% successful vs. unsuccessful object-cued location recall.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath('E:\TreasureHunt\Functions\');
addpath('E:\OpenField\Functions\');

% variables
param               = [];
param.clusterName   = 'Cluster4Analysis_LK_20200416_113912';
param.timeRes       = 0.25; % temporal resolution
param.timeLimits    = [-5, 1]; % temporal limits relative to response
param.preDur        = 1; % duration of time period before location recall used as baseline
param.postDur       = 1; % duration of time period after response
param.smoothKernel  = 'gaussian';
param.smoothFac     = 1; % smoothing factor (1 = no smoothing)
param.memoryCutoff  = 0.9; % cutoff for good memory trials
param.roundingPrec  = 3; % rounding precision
param.numSurrogates = 1001; % number of surrogates
param.myRNG         = 111;
rng(param.myRNG);

% paths
paths           = [];
paths.spike     = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.beh       = 'E:\TreasureHunt\Beh_210420\';
paths.save      = strcat('E:\TreasureHunt\LocRecall_20200617\20210418\');
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20201113\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

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

%% load previous ELRPD-cell results

% load previous ELRPD-cell results
ELRPD   = load('E:\TreasureHunt\ELRPDCellAnalysis_220420\20200722_Dir30_Loc6x6_EgoLocDir30\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', ...
    'allRes', 'param');

%% preallocations

% main results for each cell
allRes          = [];

% for sanity checking
exByWC          = []; % excluded by wave-clus
exByELRPD       = []; % excluded due to ELRPD-cell analysis

%% loop through subjects
for iSub = 1:length(subjects)
    
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
        
        %% wires to investigate
        
        % available microwires
        wires       = dir(fullfile(sessions(iSess).folder, sessions(iSess).name, 'chan*'));
        tmp         = split(transpose({wires.name}), 'n');
        [~, I]      = sort(cellfun(@str2num, tmp(:, 2)));
        wires       = wires(I);
        
        %% loop through wires
        for iWire = 1:size(wires, 1)
            
            % report
            fprintf('\tWire: %s.\n', wires(iWire).name);
            
            %% get spike times in behavioral time
            
            % load wave-clus output
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
                                
                % continue if this unit was not included in the ELRPD-cell
                % analysis
                if ~any(all([iSub, sessIdx, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2))
                    exByELRPD       = cat(1, exByELRPD, [iSub, sessIdx, iWire, iClus]); % bookkeeping
                    fprintf('This cluster was not included in the ELRPD-cell analysis.\n');
                    continue;
                end
                
                % data for this cluster only
                thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime, behavioral time
                thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
                
                %% timecourse of firing rate during location recall
                
                % onsets and ends
                onsEnds     = [trialInfo.timeCue4LocRecall, trialInfo.timeLocRecall + param.postDur];
                locRecDur   = diff(onsEnds, [], 2);
                
                % number of spikes, times, and firing rates
                respLocRecTrialTime = cell(size(onsEnds, 1), 1); % response-locked trial time
                respLocRecTrialFR   = cell(size(onsEnds, 1), 1); % response-locked FRs
                for iOns = 1:size(onsEnds, 1)
                                        
                    % response-locked timepoint-specific firing rates
                    respTrialTimeBorders            = fliplr(onsEnds(iOns, 2):-param.timeRes:onsEnds(iOns, 1)); % exactly aligned to end of recall period (plus buffer)
                    trialDurs                       = diff(respTrialTimeBorders, [], 2);
                    trialSpikes                     = histcounts(thisCluster(:, 3), respTrialTimeBorders);
                    respLocRecTrialFR{iOns, 1}      = trialSpikes ./ trialDurs;
                    trialTime                       = movmean(respTrialTimeBorders, 2, 'endpoints', 'discard');
                    respLocRecTrialTime{iOns, 1}    = trialTime - (max(respTrialTimeBorders) - param.postDur);
                end
                
                %% firing rate during a time period before location recall (baseline)
                
                % onsets and ends
                onsEnds     = [trialInfo.timeCue4LocRecall - param.preDur, trialInfo.timeCue4LocRecall];
                preRecDur   = diff(onsEnds, [], 2);
                
                % number of spikes
                numSpikes   = nan(size(onsEnds, 1), 1);
                for iOns = 1:size(onsEnds, 1)
                    
                    % total number of spikes in this period
                    numSpikes(iOns, :)  = histcounts(thisCluster(:, 3), onsEnds(iOns, :));
                end
                
                % mean firing rate
                preRecFR    = numSpikes ./ preRecDur;
                
                %% collect information for this unit
                
                % basics
                unitRes                     = [];
                unitRes.idx                 = [iSub, sessIdx, iWire, iClus];
                                
                % time-resolved response-locked firing rates
                unitRes.respLocRecTrialFR   = respLocRecTrialFR;
                unitRes.respLocRecTrialTime = respLocRecTrialTime;
                                
                % average firing rates during the baseline period
                unitRes.preRecFR            = preRecFR;
                unitRes.preRecDur           = preRecDur;
                
                % trial information
                unitRes.trialInfo           = trialInfo;
                unitRes.locRecDur           = locRecDur;
                
                % collect across units
                allRes  = cat(1, allRes, unitRes);
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load previously saved results

% load results
r   = load(strcat(paths.save, 'results.mat'));
fprintf('\nTotal number of cells: %d.\n', size(r.allRes, 1));

% indices of all units
allUnitIdx  = cell2mat({r.allRes.idx}');

% ELRPD cells
bELRPDCell = cell2mat({r.ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELRPDCell));

%% direction cells and place-like cells

% identify direction- and place-like cells
PxD             = load('E:\TreasureHunt\PlaceDirAnalysis_220520\20200623_Dir30_Loc6x6\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', 'allRes');
bDirCell      	= cell2mat({PxD.allRes.dir_Frank}') > 0.95;
bPlaceCell    	= cell2mat({PxD.allRes.loc_Frank}') > 0.95;

% cells that are either ELRPD cells or direction cells or place cells
bSpatialCell   	= bELRPDCell | bDirCell | bPlaceCell;

%% memory performance in each trial

% random locations within the arena
cfg         = [];
cfg.maxR    = ELRPD.param.arenaRadius;
cfg.minR    = 0;
cfg.N       = 1000001;
cfg.centerX = ELRPD.param.arenaCtr(1);
cfg.centerY = ELRPD.param.arenaCtr(2);
randLocs    = LK_RandomPointsInCircle_101119(cfg);

% loop through cells
allMemPerf  = cell(size(r.allRes, 1), 1);
for iCell = 1:size(r.allRes, 1)
    
    % compute memory performance in each trial if it is a new subject or a
    % new session
    if iCell == 1 || allUnitIdx(iCell, 1) ~= allUnitIdx(iCell - 1, 1) || allUnitIdx(iCell, 2) ~= allUnitIdx(iCell - 1, 2)
        
        % trial-wise drop error and correction for potential drop errors
        dropErrors      = sqrt((r.allRes(iCell).trialInfo.CHOSEN_TEST_POSITION(:, 1) - r.allRes(iCell).trialInfo.CORRECT_TEST_POSITION(:, 1)) .^ 2 + ...
            (r.allRes(iCell).trialInfo.CHOSEN_TEST_POSITION(:, 3) - r.allRes(iCell).trialInfo.CORRECT_TEST_POSITION(:, 3)) .^ 2);
        surroDropErrors = pdist2(r.allRes(iCell).trialInfo.CORRECT_TEST_POSITION(:, [1, 3]), randLocs);
        memPerf         = sum(dropErrors < surroDropErrors, 2) ./ size(surroDropErrors, 2);
    end
    
    % collect across cells
    allMemPerf{iCell, 1}    = memPerf; % memory performance in each trial, stored separately for each cell
end

%% investigate time-resolved FRs

% trials with good memory performance, stored for each cell
allBGoodMem             = cell(size(r.allRes, 1), 1);

% response-locked firing rates
respTrialTime           = round(min(cellfun(@min, lk_unwrapNestedCell({r.allRes.respLocRecTrialTime}'))):r.param.timeRes:...
    max(cellfun(@max, lk_unwrapNestedCell({r.allRes.respLocRecTrialTime}'))), param.roundingPrec); % minimum time of any trial until maximum time of any trial
respLocRecTrialFR       = nan(size(r.allRes, 1), numel(respTrialTime)); % all trials (averaged within each cell)
respLocRecTrialFR_good  = nan(size(r.allRes, 1), numel(respTrialTime)); % good trials (averaged within each cell)
respLocRecTrialFR_bad   = nan(size(r.allRes, 1), numel(respTrialTime)); % bad trials (averaged within each cell)

% control analyses
numDataPoints           = nan(size(r.allRes, 1), numel(respTrialTime)); % count for each cell and time point how many trials contributed
numTrials_goodVSbad     = nan(size(r.allRes, 1), 2); % number of good trials and bad trials

% loop through cells
for iCell = 1:size(r.allRes, 1)
    
    % determine trials with good memory performance
    bGoodMem                        = allMemPerf{iCell, 1} >= param.memoryCutoff;
    allBGoodMem{iCell}              = bGoodMem; % store for each cell
    numTrials_goodVSbad(iCell, :) 	= [sum(bGoodMem), sum(~bGoodMem)];
    
    % loop through trials to align the time-resolved FRs
    respOrigTrialFR = nan(size(r.allRes(iCell).respLocRecTrialFR, 1), numel(respTrialTime)); % response-locked
    for iTrial = 1:size(r.allRes(iCell).respLocRecTrialFR, 1)
                
        % put FRs onto common timeline, response locked
        tmpFR                           = r.allRes(iCell).respLocRecTrialFR{iTrial, 1};
        bcFR                            = tmpFR - r.allRes(iCell).preRecFR(iTrial, 1); % subtract pre-recall FR
        logIdx                          = ismember(respTrialTime, round(r.allRes(iCell).respLocRecTrialTime{iTrial}, param.roundingPrec));
        respOrigTrialFR(iTrial, logIdx) = smoothdata(bcFR, 2, param.smoothKernel, r.param.smoothFac); % smooth along second dimension
    end
    
    % original, absolute timeline, response-locked
    numDataPoints(iCell, :)             = sum(~isnan(respOrigTrialFR), 1); % sum across trials
    respLocRecTrialFR(iCell, :)         = nanmean(respOrigTrialFR, 1); % average across trials
    respLocRecTrialFR_good(iCell, :)    = nanmean(respOrigTrialFR(bGoodMem, :), 1);
    respLocRecTrialFR_bad(iCell, :)     = nanmean(respOrigTrialFR(~bGoodMem, :), 1);
end

%% fieldtrip

% add fieldtrip to path
addpath(paths.fieldtrip);
ft_defaults;
ft_warning off;

%% cluster-based permutation testing and figure for response-locked firing rates, good vs. bad

fprintf('\n=== Response-locked firing rates, good vs. bad trials.\n');

% reset rng
rng(r.param.myRNG);

% data groups and time for plotting
groups  = {'spatialCells'; 'nonSpatialCells'; 'ELRPDCells'};
TWOI    = respTrialTime >= param.timeLimits(1) & respTrialTime <= param.timeLimits(2); % time window of interest
time  	= respTrialTime(TWOI);
myXLim  = param.timeLimits; % for plotting
myYLim  = [-0.5, 0.7]; % for plotting

% loop through groups
for iGroup = 1:size(groups, 1)
    
    % specific data
    fprintf('\nGroup: %s.\n', groups{iGroup});
    if strcmp(groups{iGroup}, 'spatialCells')
        bCellMask   = bSpatialCell;
    elseif strcmp(groups{iGroup}, 'nonSpatialCells')
        bCellMask   = ~bSpatialCell;
    elseif strcmp(groups{iGroup}, 'ELRPDCells')
        bCellMask   = bELRPDCell;
    end
    dataGood    = respLocRecTrialFR_good(bCellMask, TWOI);
    dataBad     = respLocRecTrialFR_bad(bCellMask, TWOI);
    
    %% cluster-based permutation testing using fieldtrip, good > bad
    
    % re-organize data for fieldtrip
    timelock1   = cell(1, size(dataGood, 1));
    timelock2  	= cell(1, size(dataBad, 1));
    for iCell = 1:size(dataGood, 1)
        
        % good trials
        timelock1{iCell}.avg       = dataGood(iCell, :);
        timelock1{iCell}.label     = {'neuron'};
        timelock1{iCell}.time      = time;
        
        % bad trials
        timelock2{iCell}.avg       = dataBad(iCell, :);
        timelock2{iCell}.label     = {'neuron'};
        timelock2{iCell}.time      = time;
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = 'all';
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'depsamplesT';
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum';
    cfg.neighbours          = [];
    cfg.tail                = 1; % H1: greater firing rates for good > bad trials
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = param.numSurrogates;
    numCells                = size(dataGood, 1);
    design                  = zeros(2, numCells * 2);
    design(1, :)            = [1:numCells, 1:numCells]; % cell indices
    design(2, :)            = [ones(1, numCells), ones(1, numCells) * 2]; % 1 = good; 2 = bad
    cfg.design              = design; % design matrix
    cfg.uvar                = 1; % row of the design matrix that contains the units of observation
    cfg.ivar                = 2; % row of the design matrix that contains the independent variable
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    
    % report p-value of positive cluster
    if isfield(outFT, 'posclusters') && ~isempty(outFT.posclusters)
        thisProb    = cell2mat({outFT.posclusters.prob});
        fprintf('\nFieldtrip evaluation. P-value of positive cluster: p = %.6f.\n', thisProb(thisProb < 0.05));
    end
    
    %% cluster-based permutation testing, good > 0
    
    % test for significant activation above zero
    c             	= [];
    c.mat           = dataGood;
    c.alpha         = 0.05;
    c.direction     = 'pos'; % H1: firing rates are above zero
    c.numSurrogates = param.numSurrogates;
    outCBPT         = LK_PermutationTest_OneSampleAgainstZero_20200708(c);
        
    %% figure
    
    % mean and SEM
    mGood   = nanmean(dataGood);
    steGood = nanstd(dataGood) ./ sqrt(sum(~isnan(dataGood), 1));
    mBad    = nanmean(dataBad);
    steBad  = nanstd(dataBad) ./ sqrt(sum(~isnan(dataBad), 1));
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([time, fliplr(time)], ...
        [mBad - steBad, fliplr(mBad + steBad)], ...
        [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(time, mBad, ...
        'Color', [1, 0, 0], 'LineWidth', 2);
    patch([time, fliplr(time)], ...
        [mGood - steGood, fliplr(mGood + steGood)], ...
        [0, 0.5, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(time, mGood, ...
        'Color', [0, 0.5, 0], 'LineWidth', 2);
    % significance info
    LK_SigLine(time, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT.mask);
    LK_SigLine(time, [max(myYLim) - range(myYLim) * 0.04; max(myYLim) - range(myYLim) * 0.065], outCBPT.logIdxSigClus, [0.5, 0.5, 0.5]);
    % adjust axis
    set(gca, ...
        'xlim', [min(myXLim), max(myXLim)], 'xtick', min(myXLim):max(myXLim), 'ylim', myYLim, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    % adjust font
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(r.paths.save, groups{iGroup}, '_FRsOnOriginalTime_respLocked_goodVSbadMemory_20201115'), '-dtiff', '-r450');
end

%% cluster-based permutation test for interaction (cell type x performance) using fieldtrip, response-locked time

fprintf('\n=== Testing for an interaction between cell type and performance, using cluster-based permutation testing across the entire response-locked time window.\n');

% groups of cells to contrast
groups  = {'spatialCells', 'nonSpatialCells'; ...
    'ELRPDCells', 'nonSpatialCells'};

% time window of interest
TWOI        = respTrialTime >= param.timeLimits(1) & respTrialTime <= param.timeLimits(2);
thisTime    = respTrialTime(TWOI);

% loop through groups
for iGroup = 1:size(groups, 1)
    
    % report
    fprintf('\n... Interaction good vs bad trials and %s vs %s.\n\n', groups{iGroup, 1}, groups{iGroup, 2});
    
    % difference between conditions, separately for the two cell types
    if strcmp(groups{iGroup, 1}, 'ELRPDCells')
        groupA_goodVSbad 	= respLocRecTrialFR_good(bELRPDCell, TWOI) - respLocRecTrialFR_bad(bELRPDCell, TWOI); % ELRPD cells, good - bad
    elseif strcmp(groups{iGroup, 1}, 'spatialCells')
        groupA_goodVSbad    = respLocRecTrialFR_good(bSpatialCell, TWOI) - respLocRecTrialFR_bad(bSpatialCell, TWOI); % spatial cells, good - bad
    end
    if strcmp(groups{iGroup, 2}, 'nonSpatialCells')
        groupB_goodVSbad    = respLocRecTrialFR_good(~bSpatialCell, TWOI) - respLocRecTrialFR_bad(~bSpatialCell, TWOI); % non-spatial cells, good - bad
    end
    
    % re-organize data for fieldtrip
    timelock1       = [];
    timelock1.label = {'neuron'};
    for iCell = 1:size(groupA_goodVSbad, 1)
        timelock1.time{1, iCell}    = thisTime;
        timelock1.trial{1, iCell}   = groupA_goodVSbad(iCell, :);
    end
    timelock2       = [];
    timelock2.label = {'neuron'};
    for iCell = 1:size(groupB_goodVSbad, 1)
        timelock2.time{1, iCell}    = thisTime;
        timelock2.trial{1, iCell}   = groupB_goodVSbad(iCell, :);
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = 'all';
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'indepsamplesT';
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum';
    cfg.neighbours          = [];
    cfg.tail                = 1; % H1: greater firing rates for spatial/ELRPD cells as compared to non-spatial cells
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = param.numSurrogates;
    cfg.design              = [ones(size(timelock1.trial)), 2 .* ones(size(timelock2.trial))]; % design matrix
    cfg.ivar                = 1;
    
    % fieldtrip estimation
    outFT_inter = ft_timelockstatistics(cfg, timelock1, timelock2);
    
    % report p-value of largest positive cluster
    fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT_inter.posclusters.prob})));
    fprintf('Time window of the largest positive cluster: %.3f to %.3f.\n', min(thisTime(outFT_inter.mask)), max(thisTime(outFT_inter.mask)));
        
    %% figure for interaction effect
    
    % mean and standard error
    m   = {nanmean(groupB_goodVSbad, 1), nanmean(groupA_goodVSbad, 1)}; % difference #2, difference #1
    sem = {nanstd(groupB_goodVSbad, 1) ./ sqrt(sum(~isnan(groupB_goodVSbad), 1)), ...
        nanstd(groupA_goodVSbad, 1) ./ sqrt(sum(~isnan(groupA_goodVSbad), 1))}; % difference #2, difference #1
    
    % colors
    myColor1    = [0.7, 0.7, 0.7];
    myColor2    = [0.4, 0.4, 0.4];
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([thisTime, fliplr(thisTime)], ...
        [m{1} - sem{1}, fliplr(m{1} + sem{1})], myColor1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    patch([thisTime, fliplr(thisTime)], ...
        [m{2} - sem{2}, fliplr(m{2} + sem{2})], myColor2, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(thisTime, m{1}, '-', ...
        'Color', myColor1, 'LineWidth', 2);
    plot(thisTime, m{2}, '-', ...
        'Color', myColor2, 'LineWidth', 2);
    set(gca, ...
        'YLim', myYLim);
    % significance info
    LK_SigLine(thisTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT_inter.mask);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(myXLim), max(myXLim)], 'xtick', min(myXLim):max(myXLim), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(paths.save, 'LocRecall_InteractionEffectCBPT_goodVSbad_', groups{iGroup, 1}, 'VS', groups{iGroup, 2}), '-dtiff', '-r600');  
end
