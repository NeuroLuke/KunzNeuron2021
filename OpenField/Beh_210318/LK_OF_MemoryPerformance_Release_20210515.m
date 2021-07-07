%==========================================================================
% Analysis of spatial memory performance.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clc; clear; close all;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param           = [];
param.myRNG     = 111;
rng(param.myRNG);

% paths
paths           = [];
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.save      = paths.beh;

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

% possible random locations within arena
dt              = [];
dt.maxR       	= 5000; % maximum radius
dt.minR     	= 0; % minimum radius
dt.N          	= 1000000; % number of locations to create
dt.centerX    	= 0; % arena center x
dt.centerY    	= 0; % arena center y
locsInArena   	= LK_RandomPointsInCircle_101119(dt);

% preallocate memory accuracy
allMemAcc                   = cell(size(subjects, 1), 1);
allMemAcc_trialFirstLast    = nan(size(subjects, 1), 2);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    %% load data
    
    % load trial information
    behFile  	= dir(strcat(paths.beh, subjects{iSub, 1}, '\trials.mat'));
    tmp       	= load(strcat(behFile.folder, filesep, behFile.name));
    trials    	= tmp.trials;
    
    % exclude trials that were not fully completed
    trials    	= trials(~isnan(trials(:, 8)), :);
    
    %% drop error and memory accuracy/performance
    
    % drop error and memory accuracy per trial
    dropError   = trials(:, end); % drop error
    memAcc      = nan(size(dropError)); % memory accuracy
    for iTrial = 1:size(memAcc, 1)
        
        % object location in this trial
        objLoc          = trials(iTrial, 9:10);
        
        % distance of this object location to other arena locations
        D2locsInArena   = pdist2(objLoc, locsInArena);
        
        % estimate memory accuracy (1 = perfect)
        memAcc(iTrial, 1)   = sum(dropError(iTrial) < D2locsInArena) / numel(D2locsInArena);
    end
    
    % collect across subjects
    allMemAcc{iSub, 1}      = memAcc;
    
    %% memory performance in first and last trial
    
    % memory performance in first and last trial, averaged across the eight
    % different objects
    objects                 = trials(:, 2);
    memAcc_trialFirstLast   = nan(numel(unique(objects)), 2);
    for iObj = 1:numel(unique(objects))
        thisObjMemAcc                   = memAcc(objects == (iObj - 1)); % object indices range from 0:7
        memAcc_trialFirstLast(iObj, :)  = [thisObjMemAcc(1), thisObjMemAcc(end)];
    end
    allMemAcc_trialFirstLast(iSub, :)   = mean(memAcc_trialFirstLast, 1); % average across different objects    
end

%% results

% percent of trials with performance above chance
fprintf('Percentage of trials with a normalized memory accuracy above 0.5: %.3f%%.\n', 100 * sum(cell2mat(allMemAcc) > 0.5) / numel(cell2mat(allMemAcc)));

%% overall histogram showing memory accuracy per trial

% histogram of all memory accuracies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
hold on;
myH = histogram(cell2mat(allMemAcc), 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
plot([0.5, 0.5], [0, max(myH.Values)], ':', ...
    'Color', [1, 0, 0], 'LineWidth', 2);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [-0.05, 1.05], ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Memory performance');
yl = ylabel('# trials', ...
    'units', 'normalized', 'position', [-0.05, 0.5, 0]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'AllTrials_memAcc_020220'), '-dtiff', '-r450');

%% memory performance in the last vs. first trial

% compare memory accuracy in the first vs. last trial
[~, p, ~, stats]   = ttest(allMemAcc_trialFirstLast(:, end), allMemAcc_trialFirstLast(:, 1));
fprintf('Testing memory accuracy in the last vs. first trial: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.55, 3.8, 4]);
hold on;
for iSub = 1:size(subjects, 1)
    plot([1, 2], allMemAcc_trialFirstLast(iSub, :), '-', ...
        'Color', rgb('gray'));
end
plot([1, 2], mean(allMemAcc_trialFirstLast), '-', ...
    'Color', rgb('blue'), 'LineWidth', 3); % average across sessions
hold off;
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'First', 'Last'}, ...
    'ylim', [0, 1], 'ytick', [0, 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
xlabel('Trial', ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
ylabel('Memory performance', ...
    'units', 'normalized', 'position', [-0.15, 0.5, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'AllSubj_memAcc_TrialFirstLast_020220'), '-dtiff', '-r450');

%% memory performance in the last vs. first trial, separetely for first and second sessions

% identify first-session subjects
minLength   = min(cellfun(@(x) numel(x), subjects));
cutSubjects = cell(size(subjects, 1), 1);
for iSub = 1:size(subjects, 1)
    cutSubjects{iSub, 1}    = subjects{iSub}(1:minLength);
end
[~, uniqueIdx]  = unique(cutSubjects);

% subjects from the first vs. second session
b1stSess        = ismember(subjects, subjects(uniqueIdx));

% compare memory accuracy in first vs. last trial, first sessions
[~, p, ~, stats]    = ttest(allMemAcc_trialFirstLast(b1stSess, end), allMemAcc_trialFirstLast(b1stSess, 1));
fprintf('Testing memory accuracy in last vs. first trial, 1st sessions only: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% compare memory accuracy in first vs. last trial, second sessions
[~, p, ~, stats]    = ttest(allMemAcc_trialFirstLast(~b1stSess, end), allMemAcc_trialFirstLast(~b1stSess, 1));
fprintf('Testing memory accuracy in last vs. first trial, 2nd sessions only: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
