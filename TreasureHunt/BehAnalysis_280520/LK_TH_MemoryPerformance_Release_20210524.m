%==========================================================================
% This script examines the subjects' memory performance in TreasureHunt.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clc; close all; clear;
addpath('E:\TreasureHunt\Functions\');
addpath('E:\OpenField\Functions\');
rng(1);

% paths
paths               = [];
paths.beh           = 'E:\TreasureHunt\Beh_210420\';
paths.spike         = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.save          = 'E:\TreasureHunt\BehAnalysis_280520\';

% settings
param             	= [];
param.arenaXLim   	= [430, 310]; % left-to-right
param.arenaZLim    	= [420, 300]; % bottom-to-top
param.arenaCtr     	= [mean(param.arenaXLim), mean(param.arenaZLim)]; % arena center (x/z)
param.arenaRadius  	= 50; % arena radius

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

% random locations within the arena
cfg         = [];
cfg.maxR    = param.arenaRadius;
cfg.minR    = 0;
cfg.N       = 10000000;
cfg.centerX = param.arenaCtr(1);
cfg.centerY = param.arenaCtr(2);
randLocs    = LK_RandomPointsInCircle_101119(cfg);

% preallocate
objRecallPerf   = cell(size(subjects, 1), 1);
locRecallPerf   = cell(size(subjects, 1), 1);
objRecallDur    = cell(size(subjects, 1), 1);
locRecallDur    = cell(size(subjects, 1), 1);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % available sessions
    sessions    = dir(strcat(paths.spike, subjects{iSub}, '\session_*'));
    
    % preallocate performance for this subject
    subjObjRecallPerf   = nan(size(sessions, 1), 1);
    subjLocRecallPerf   = cell(size(sessions, 1), 1);
    subjObjRecallDur    = cell(size(sessions, 1), 1);
    subjLocRecallDur    = cell(size(sessions, 1), 1);
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % report progress
        fprintf('Working on subject "%s", session "%s".\n', subjects{iSub}, sessions(iSess).name);
        
        % trial information for this session
        tmp         = load(fullfile(paths.beh, subjects{iSub}, filesep, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo   = tmp.trialInfo;
        
        % memory performance during object recall
        subjObjRecallPerf(iSess, 1)     = sum(strcmp(trialInfo.CORTANA_RESPONSE, 'CORRECT')) / size(trialInfo.CORTANA_RESPONSE, 1);
        
        % memory performance during location recall
        thisLocRecallPerf   = nan(size(trialInfo.CORRECT_TEST_POSITION, 1), 1);
        for iTrial = 1:size(trialInfo.CORRECT_TEST_POSITION, 1)
            dropError                       = pdist2(trialInfo.CHOSEN_TEST_POSITION(iTrial, [1, 3]), trialInfo.CORRECT_TEST_POSITION(iTrial, [1, 3]));
            potDropErrors                   = pdist2(trialInfo.CORRECT_TEST_POSITION(iTrial, [1, 3]), randLocs);
            thisLocRecallPerf(iTrial, 1)    = sum(dropError < potDropErrors) / numel(potDropErrors);
        end
        subjLocRecallPerf{iSess, 1} = thisLocRecallPerf;
        
        % duration of object recall
        subjObjRecallDur{iSess, 1}  = trialInfo.RECORDING_ENDED - trialInfo.RECORDING_STARTED;
        
        % duration of location recall
        subjLocRecallDur{iSess, 1}  = trialInfo.timeLocRecall - trialInfo.timeCue4LocRecall;
    end
    
    % collapse across subjects
    objRecallPerf{iSub, 1}  = subjObjRecallPerf; % session-wise average values
    locRecallPerf{iSub, 1}  = subjLocRecallPerf; % trial-wise values
    objRecallDur{iSub, 1}   = subjObjRecallDur;
    locRecallDur{iSub, 1}   = subjLocRecallDur;
end

%% figure for performance during (location-cued) object recall

% histogram of subject-wise memory accuracies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 5]);
axes('units', 'centimeters', 'position', [1.5, 2.2, 3.8, 2.5]);
hold on;
myH = histogram(cell2mat(objRecallPerf), 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [-0.05, 1.05], ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel({'Object recall', 'memory performance'}, ...
    'units', 'normalized', 'position', [0.5, -0.34, 0]);
yl = ylabel('# sessions', ...
    'units', 'normalized', 'position', [-0.2, 0.5, 0]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'AllSubj_ObjectRecall_MemoryPerformance'), '-dtiff', '-r450');

%% figure for performance during (object-cued) location recall

% trial-wise location-recall performance
trialLocRecallPerf  = cell2mat(lk_unwrapNestedCell(locRecallPerf));

% histogram of trial-wise memory accuracies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 5]);
axes('units', 'centimeters', 'position', [1.5, 2.2, 3.8, 2.5]);
hold on;
myH = histogram(trialLocRecallPerf, 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
plot([0.5, 0.5], [0, max(myH.Values)], ':', ...
    'Color', [1 0 0], 'LineWidth', 2);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [-0.05, 1.05], ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel({'Location recall', 'memory performance'}, ...
    'units', 'normalized', 'position', [0.5, -0.34, 0]);
yl = ylabel('# trials', ...
    'units', 'normalized', 'position', [-0.2, 0.5, 0]);
set(gca, ...
    'ylim', [0, max(myH.Values)], 'ytick', [0, max(myH.Values)]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'AllSubj_LocationRecall_MemoryPerformance'), '-dtiff', '-r450');

%% duration of location recall

% durations of all location recalls
allLocRecallDur         = cell2mat(lk_unwrapNestedCell(locRecallDur));

% average duration of location recalls per session
sessLocRecallDur        = cellfun(@mean, lk_unwrapNestedCell(locRecallDur));

% report
meanSessLocRecallDur    = mean(sessLocRecallDur);
semSessLocRecallDur     = std(sessLocRecallDur) ./ sqrt(size(sessLocRecallDur, 1));
fprintf('Average duration of location recall: %.3f +/- %.3f sec (grand min = %.3f sec; grand max = %.3f sec).\n', meanSessLocRecallDur, semSessLocRecallDur, ...
    min(allLocRecallDur), max(allLocRecallDur));
