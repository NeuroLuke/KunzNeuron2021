%==========================================================================
% This script examines the timecourse of object cells and egocentric
% bearing (ELRPD) cells during spatial memory recall during cue
% presentation.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));
addpath('E:\TreasureHunt\Functions\');

% settings
param                   = [];
param.myRNG             = 555;
param.clusterName       = 'HardCrit_20190618_091853';
param.timeRes           = 0.25; % temporal resolution
param.smoothFac         = 1; % smoothing factor (1 = no smoothing)
param.trialTimeBorders  = -1:param.timeRes:3; % relative to cue onset
param.trialTime         = movmean(param.trialTimeBorders, 2, 'endpoints', 'discard');
param.baseTimeBorders   = [-1, 0]; % baseline time period
param.cutoffQuantile    = 2/8; % cutoff for differentiating between close and far objects

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20201113\';
paths.save      = strcat('E:\OpenField\CuePresentation_20200803\20210515\');
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

% load reference analysis (to only include egocentric bearing cells in the
% analysis)
ELRPD   = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'locDir');

%% preallocations

% main results for each cell
allRes     	= [];

% bookkeeping
exByWC    	= []; % excluded by wave-clus
exByRefAna 	= []; % excluded by reference analysis

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(param.myRNG);
    
    %% behavioral data
    
    % load "trials" for relevant information
    tmp     = load(strcat(paths.beh, subjects{iSub}, filesep, 'trials.mat')); % trialidx, objectidx, ITI, cue, retrieval, feedback, reencoding, grab, correctXY, dropXY, drop error
    trials  = tmp.trials;
    
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
        
        %% get spike times in behavioral time
        
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
            error('Size of "cluster_class" not congruent with "cluster_class_behtime');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
            
            % continue if this cluster was not included in the
            % egocentric-bearing-cell analysis
            refIdx  = all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2);
            if sum(refIdx) ~= 1
                fprintf('\t\t- This cluster was not included in the reference analysis.\n');
                exByRefAna  = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster     = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            
            %% firing rate during cue presentation
            
            % trial onsets
            trialOns            = trials(:, 4);
            trialTimeBorders    = trialOns + param.trialTimeBorders;
            durations           = diff(trialTimeBorders, [], 2);
            
            % time-resolved FRs during cue
            numSpikes   = nan(size(durations));
            for iOns = 1:size(trialOns, 1)
                numSpikes(iOns, :)  = histcounts(thisCluster(:, 3), trialTimeBorders(iOns, :));
            end
            
            % time-resolved FRs during cue
            cueFR   = numSpikes ./ durations;
            
            %% collect information for this unit
            
            % basics
            unitRes      	= [];
            unitRes.idx   	= [iSub, iWire, iClus];
            
            % time-resolved FRs
            unitRes.cueFR 	= cueFR;
            
            % behavior
            unitRes.trials  = trials;
            
            % collect across units
            allRes          = cat(1, allRes, unitRes);
        end
    end
end

%% save all results

% save
save(strcat(paths.save, 'results'));

%% load previously saved results

% load results
r   = load(strcat(paths.save, 'results.mat'));
fprintf('Total number of cells: %d.\n', size(r.allRes, 1));

% indices of all units
allUnitIdx  = cell2mat({r.allRes.idx}');

%% cell types

% egocentric bearing cells
bELRPDCell      = cell2mat({ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of egocentric bearing cells: %d.\n', sum(bELRPDCell));

% object cells
obj             = load('E:\OpenField\ObjectCellAnalysis_010120\20200803_Dir30_Loc10x10_Obj\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat');
allObjFrank     = cell2mat({obj.allRes.obj_Frank}');
numPrefObj      = sum(cell2mat({obj.allRes.objBins})', 2);
bObjCell        = allObjFrank > 0.95 & numPrefObj >= 1;
fprintf('Number of object cells: %d.\n', sum(bObjCell));

%% time-resolved firing rates during cue presentation

% extract average time-resolved firing rates
allMeanCueFR_prefObj    = nan(size(r.allRes, 1), numel(r.param.trialTime)); % preferred objects
allMeanCueFR_unprefObj  = nan(size(r.allRes, 1), numel(r.param.trialTime)); % unpreferred objects
allMeanCueFR_close2COM  = nan(size(r.allRes, 1), numel(r.param.trialTime)); % close to reference point
allMeanCueFR_far2COM    = nan(size(r.allRes, 1), numel(r.param.trialTime)); % far from reference point
for iCell = 1:size(r.allRes, 1)

    % smoothing and baseline firing rate
    thisFR      = r.allRes(iCell).cueFR;
    thisFR      = smoothdata(thisFR, 2, 'gaussian', r.param.smoothFac); % smooth along 2nd dimension
    baselineFR  = mean(thisFR(:, r.param.trialTime > param.baseTimeBorders(1) & r.param.trialTime < param.baseTimeBorders(2)), 2);
    
    % exclude first and last trial (because, due to preprocessing, there
    % are no spikes before the onset of the first cue and after the onset
    % of the last cue)
    bValid      = r.allRes(iCell).trials(:, 1) ~= 1 & r.allRes(iCell).trials(:, 1) ~= max(r.allRes(iCell).trials(:, 1));
    
    % mean across trials with preferred object
    prefObjTrials                   = ismember(r.allRes(iCell).trials(:, 2), obj.object.idx(obj.allRes(iCell).objBins == 1));
    allMeanCueFR_prefObj(iCell, :)  = mean(thisFR(prefObjTrials & bValid, :) - baselineFR(prefObjTrials & bValid, :));
    
    % mean across trials with unpreferred object
    allMeanCueFR_unprefObj(iCell, :)    = mean(thisFR(~prefObjTrials & bValid, :) - baselineFR(~prefObjTrials & bValid, :));
    
    % differentiate between trials in which the object location is close
    % vs. far from the reference point
    thisCOMxy        	= ELRPD.allRes(iCell).locdir_COMxy;
    thisRespLocs        = r.allRes(iCell).trials(:, 9:10);
    D_respLocs2COMxy    = pdist2(thisRespLocs, thisCOMxy);
    cutoff              = quantile(unique(D_respLocs2COMxy(bValid)), r.param.cutoffQuantile);
    bRespClose2COMxy    = D_respLocs2COMxy < cutoff;
        
    % mean across trials with response location close to reference point
    allMeanCueFR_close2COM(iCell, :)    = mean(thisFR(bRespClose2COMxy & bValid, :) - baselineFR(bRespClose2COMxy & bValid, :));
    
    % mean across trials with response location far from reference point
    allMeanCueFR_far2COM(iCell, :)      = mean(thisFR(~bRespClose2COMxy & bValid, :) - baselineFR(~bRespClose2COMxy & bValid, :));
end

%% evaluation: object cells during cue presentation, pref. vs. unpref. objects (paired test)

% report
fprintf('\n==== Testing firing rates during pref. vs. unpref. objects, object cells.\n\n');

% add fieldtrip to path
addpath(paths.fieldtrip);
ft_defaults;
ft_warning off;

% select data
selDataPref     = allMeanCueFR_prefObj(bObjCell, :);
selDataUnpref   = allMeanCueFR_unprefObj(bObjCell, :);

% reformat for fieldtrip
timelock1   = cell(1, size(selDataPref, 1));
timelock2  	= cell(1, size(selDataUnpref, 1));
for iCell = 1:size(selDataPref, 1)
        
    % condition #1
    timelock1{iCell}.avg       = selDataPref(iCell, :);
    timelock1{iCell}.label     = {'neuron'};
    timelock1{iCell}.time      = r.param.trialTime;
    
    % condition #2
    timelock2{iCell}.avg       = selDataUnpref(iCell, :);
    timelock2{iCell}.label     = {'neuron'};
    timelock2{iCell}.time      = r.param.trialTime;
end

% fieldtrip configuration
cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT'; % type of test at the first level
cfg.correctm            = 'cluster'; % correction method
cfg.clusteralpha        = 0.05; % first-level alpha used for thresholding
cfg.clusterstatistic    = 'maxsum'; % type of cluster test statistic
cfg.neighbours          = [];
cfg.tail                = 1; % H1: stronger activation during preferred vs. unpreferred trials
cfg.alpha               = 0.05; % alpha level for the permutation test
cfg.correcttail         = 'alpha';
cfg.numrandomization    = 10001; % number of surrogates
numCells                = size(selDataPref, 1);
design                  = zeros(2, numCells * 2);
design(1, :)            = [1:numCells, 1:numCells]; % cell indices
design(2, :)            = [ones(1, numCells), ones(1, numCells) * 2]; % group indices
cfg.design              = design; % design matrix
cfg.uvar                = 1; % row of the design matrix that contains the units of observation
cfg.ivar                = 2; % row of the design matrix that contains the independent variable

% fieldtrip estimation
outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

% report p-value of largest positive cluster
if isfield(outFT, 'posclusters')
    fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT.posclusters.prob})));
end

%% figure: object cells during cue presentation, pref. vs. unpref. objects (paired test)

% mean and standard error
m   = {nanmean(selDataPref, 1), nanmean(selDataUnpref, 1)}; % mean
sem = {nanstd(selDataPref, [], 1) ./ sqrt(sum(~isnan(selDataPref), 1)), ...
    nanstd(selDataUnpref, [], 1) ./ sqrt(sum(~isnan(selDataUnpref), 1))}; % SEM

% create figure
myYLim  = [-1, 3.5];
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
    [m{1} - sem{1}, fliplr(m{1} + sem{1})], rgb('orange'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3); % preferred
patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
    [m{2} - sem{2}, fliplr(m{2} + sem{2})], [0.7, 0.7, 0.7], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3); % unpreferred
plot(r.param.trialTime, m{1}, '-', ...
    'Color', rgb('orange'), 'LineWidth', 2); % preferred
plot(r.param.trialTime, m{2}, '-', ...
    'Color', [0.7, 0.7, 0.7], 'LineWidth', 2); % unpreferred
set(gca, ...
    'YLim', myYLim);
% significance info
LK_SigLine(r.param.trialTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT.mask .* outFT.posclusterslabelmat == 1);
% x- and y-axis
myxline(0, '--', [0, 0, 0], 'bottom');
myyline(0, '--', [0, 0, 0], 'bottom');
hold off;
set(gca, ...
    'xlim', [min(r.param.trialTime), max(r.param.trialTime)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Time (s)');
yl = ylabel('Firing rate (Hz)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'objCells_cueFR_prefVSUnprefObj_20201115'), '-dtiff', '-r450');

%% figure: object-egocentric-bearing cells vs. object-non-egocentric-bearing cells during cue presentation, pref. objects, unpref. objects (unpaired tests)

% report
fprintf('\n==== Testing firing rates of object-egocentric-bearing cells vs. object-non-egocentric-bearing cells, separately for pref. and unpref. objects.\n\n');

% experimental groups
groups  = {'prefObj', 'unprefObj'};
for iGroup = 1:numel(groups)
    
    % select data
    if strcmp(groups{iGroup}, 'prefObj')
        selData1    = allMeanCueFR_prefObj(bObjCell & bELRPDCell, :);
        selData2    = allMeanCueFR_prefObj(bObjCell & ~bELRPDCell, :);
        myColor     = rgb('orange');
    elseif strcmp(groups{iGroup}, 'unprefObj')
        selData1    = allMeanCueFR_unprefObj(bObjCell & bELRPDCell, :);
        selData2    = allMeanCueFR_unprefObj(bObjCell & ~bELRPDCell, :);
        myColor     = [0.7, 0.7, 0.7];
    end
    
    %% cluster-based permutation testing using fieldtrip
    
    % re-organize data for fieldtrip
    timelock1       = [];
    timelock1.label = {'neuron'};
    for iCell = 1:size(selData1, 1)
        timelock1.time{1, iCell}    = r.param.trialTime;
        timelock1.trial{1, iCell}   = selData1(iCell, :);
    end    
    timelock2       = [];
    timelock2.label = {'neuron'};
    for iCell = 1:size(selData2, 1)
        timelock2.time{1, iCell}    = r.param.trialTime;
        timelock2.trial{1, iCell}   = selData2(iCell, :);
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
    cfg.tail                = 1; % H1: stronger activation for object-egocentric-bearing cells
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = 10001;
    cfg.design              = [ones(size(timelock1.trial)), 2 .* ones(size(timelock2.trial))]; % group indices
    cfg.ivar                = 1; % row of the design matrix that contains the independent variable
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1, timelock2);
    
    % report p-value of largest positive cluster
    if isfield(outFT, 'posclusters') && ~isempty(outFT.posclusters)
        fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT.posclusters.prob})));
    end
    
    %% figure
    
    % mean and standard error
    m   = {nanmean(selData1, 1), mean(selData2, 1)}; % mean
    sem = {nanstd(selData1, [], 1) ./ sqrt(sum(~isnan(selData1), 1)), ...
        nanstd(selData2, [], 1) ./ sqrt(sum(~isnan(selData2), 1))}; % SEM
        
    % create figure
    myYLim  = [-1, 3.5];
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{1} - sem{1}, fliplr(m{1} + sem{1})], myColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
        [m{2} - sem{2}, fliplr(m{2} + sem{2})], myColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(r.param.trialTime, m{1}, '-', ...
        'Color', myColor, 'LineWidth', 2);
    plot(r.param.trialTime, m{2}, ':', ...
        'Color', myColor, 'LineWidth', 2);
    set(gca, ...
        'YLim', myYLim);
    % significance info
    LK_SigLine(r.param.trialTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT.mask .* outFT.posclusterslabelmat == 1);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(r.param.trialTime), max(r.param.trialTime)], ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(paths.save, 'objELRPDCellsVSobjNonELRPDCells_cueFR_', groups{iGroup}, '_20201115'), '-dtiff', '-r450');
end

%% figure: egocentric bearing cells during cue presentation, contrasting trials with objects close to the reference point vs. far from the reference point (paired test)

% report
fprintf('\n==== Testing firing rates of objects close to the reference point vs. far from the reference point, egocentric bearing cells.\n\n');

% select data
selDataClose    = allMeanCueFR_close2COM(bELRPDCell, :);
selDataFar      = allMeanCueFR_far2COM(bELRPDCell, :);

% re-organize data for fieldtrip
timelock1   = cell(1, size(selDataClose, 1));
timelock2  	= cell(1, size(selDataFar, 1));
for iCell = 1:size(selDataClose, 1)
    
    % close objects
    timelock1{iCell}.avg       = selDataClose(iCell, :);
    timelock1{iCell}.label     = {'neuron'};
    timelock1{iCell}.time      = r.param.trialTime;
    
    % far objects
    timelock2{iCell}.avg       = selDataFar(iCell, :);
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
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.neighbours          = [];
cfg.tail                = 1; % H1: stronger activation for close objects
cfg.alpha               = 0.05;
cfg.correcttail         = 'alpha';
cfg.numrandomization    = 10001;
numCells                = size(selDataClose, 1);
design                  = zeros(2, numCells * 2);
design(1, :)            = [1:numCells, 1:numCells]; % cell indices
design(2, :)            = [ones(1, numCells), ones(1, numCells) * 2]; % group indices
cfg.design              = design;
cfg.uvar                = 1; % row of the design matrix containing the units of observation
cfg.ivar                = 2; % row of the design matrix containing the independent variable

% fieldtrip estimation
outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

% report p-value of largest positive cluster
if isfield(outFT, 'posclusters')
    fprintf('\nFieldtrip evaluation. P-value of largest positive cluster: p = %.3f.\n', min(cell2mat({outFT.posclusters.prob})));
end

% mean and standard error
m   = {nanmean(selDataClose, 1), nanmean(selDataFar, 1)};
sem = {nanstd(selDataClose, [], 1) ./ sqrt(sum(~isnan(selDataClose), 1)), ...
    nanstd(selDataFar, [], 1) ./ sqrt(sum(~isnan(selDataFar), 1))};

% create figure
myYLim  = [-0.6, 0.8];
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
    [m{1} - sem{1}, fliplr(m{1} + sem{1})], rgb('darkgreen'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
patch([r.param.trialTime, fliplr(r.param.trialTime)], ...
    [m{2} - sem{2}, fliplr(m{2} + sem{2})], rgb('lightgreen'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(r.param.trialTime, m{1}, '-', ...
    'Color', rgb('darkgreen'), 'LineWidth', 2);
plot(r.param.trialTime, m{2}, '-', ...
    'Color', rgb('lightgreen'), 'LineWidth', 2);
set(gca, ...
    'YLim', myYLim, 'YTick', min(myYLim):0.2:max(myYLim));
if isfield(outFT, 'posclusterslabelmat')
    LK_SigLine(r.param.trialTime, [max(myYLim); max(myYLim) - range(myYLim) * 0.025], outFT.mask .* outFT.posclusterslabelmat == 1);
end
% add x- and y-axis
myxline(0, '--', [0, 0, 0], 'bottom');
myyline(0, '--', [0, 0, 0], 'bottom');
hold off;
set(gca, ...
    'xlim', [min(r.param.trialTime), max(r.param.trialTime)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Time (s)');
yl = ylabel('Firing rate (Hz)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ELRPDCells_cueFR_closeVSFar2ReferencePoint_20201115'), '-dtiff', '-r450');

%% example illustration of object locations close vs. far from the reference point

% example cell and its reference point
exCellIdx   = 679; % choose any egocentric bearing cell
thisRP      = ELRPD.allRes(exCellIdx).locdir_COMxy; % reference point

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
ax1 = axes('units', 'centimeters', 'position', [0.6, 0.6, 6.5, 6.5], 'visible', 'off');
hold on;
% object locations and their distances to the reference point
thisObjLocs     = unique(r.allRes(exCellIdx).trials(~isnan(r.allRes(exCellIdx).trials(:, 9)), 9:10), 'rows', 'stable');
D_objLocs2RP    = pdist2(thisObjLocs, thisRP);
bClose          = D_objLocs2RP <= quantile(D_objLocs2RP, r.param.cutoffQuantile);
% plot close and far objects
for iObj = 1:size(thisObjLocs, 1)
    if bClose(iObj) == true
        myColor = rgb('darkgreen');
    else
        myColor = rgb('lightgreen');
    end
    plot([thisObjLocs(iObj, 1), thisRP(1)], [thisObjLocs(iObj, 2), thisRP(2)], ':', ...
        'Color', myColor, 'LineWidth', 2);    
end
plot(thisObjLocs(~bClose, 1), thisObjLocs(~bClose, 2), 'o', ...
    'MarkerSize', 10, 'Color', rgb('lightgreen'), 'MarkerFaceColor', rgb('lightgreen'));
plot(thisObjLocs(bClose, 1), thisObjLocs(bClose, 2), 'o', ...
    'MarkerSize', 10, 'Color', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'));
% reference point
plot(thisRP(1), thisRP(2), 'o', ...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 10, 'LineWidth', 2);
% boundary
tmpX = cos(0:0.001:2*pi);
tmpY = sin(0:0.001:2*pi);
patch([tmpX .* 5000, fliplr(tmpX) .* 6250], [tmpY .* 5000, fliplr(tmpY) .* 6250], [1, 1, 1], ...
    'EdgeColor', 'none');
plot(5000 .* tmpX, 5000 .* tmpY, 'k-', ...
    'LineWidth', 2);
set(gca, ...
    'xlim', [-5400, 5400], 'ylim', [-5400, 5400], ...
    'ydir', 'reverse');
xl = text(0.5, -0.025, 'x', ...
    'units', 'normalized', 'HorizontalAlignment', 'center');
yl = text(-0.025, 0.5, 'y', ...
    'units', 'normalized', 'HorizontalAlignment', 'center');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
axis square;
% save figure
set(gcf, 'InvertHardCopy', 'off', 'Color', [1, 1, 1]);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'Illustration_objectLocsCloseVSFar2ReferencePoint'), '-dtiff', '-r450');
