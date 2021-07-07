%==========================================================================
% This script examines basic behavioral characteristics in TreasureHunt.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clc; close all; clear;
addpath('E:\TreasureHunt\Functions\');
addpath('E:\OpenField\Functions\');
rng(1);

% paths
paths       = [];
paths.beh   = 'E:\TreasureHunt\Beh_210420\';
paths.spike = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.save  = 'E:\TreasureHunt\BehAnalysis_280520\';

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

% preallocate
numSessions     = nan(size(subjects, 1), 1);
numTrials       = cell(size(subjects, 1), 1);
expDurations    = cell(size(subjects, 1), 1);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % sessions
    sessions                = dir(strcat(paths.spike, subjects{iSub}, '\session_*'));
    numSessions(iSub, 1)    = size(sessions, 1);
    
    % preallocate
    thisNumTrials       = nan(size(sessions, 1), 1);
    thisExpDurations    = nan(size(sessions, 1), 1);
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % report progress
        fprintf('Working on subject "%s", session "%s".\n', subjects{iSub}, sessions(iSess).name);
        
        % trial information for this session
        tmp         = load(fullfile(paths.beh, subjects{iSub}, filesep, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo   = tmp.trialInfo;
        
        % behavioral info
        tmp         = load(fullfile(paths.beh, subjects{iSub}, filesep, sessions(iSess).name, 'behInfo.mat'));
        behInfo     = tmp.behInfo;
        
        % number of trials
        thisNumTrials(iSess, 1)     = size(trialInfo.trialIdx, 1);
        
        % experiment duration
        thisExpDurations(iSess, 1)  = max(cell2mat(behInfo(:, 1))) - min(cell2mat(behInfo(:, 1)));
    end
    
    % collect across subjects
    numTrials{iSub, 1}      = thisNumTrials;
    expDurations{iSub, 1}   = thisExpDurations;
end

%% report

% number of sessions
fprintf('\nTotal number of sessions: %d.\n', sum(numSessions));

% number of trials/session
fprintf('Number of trials per session: %d to %d.\n', min(cell2mat(numTrials)), max(cell2mat(numTrials)));

% experiment duration
fprintf('Experiment duration per session: %.3f to %.3f minutes.\n', min(cell2mat(expDurations)) / 60, max(cell2mat(expDurations)) / 60);
