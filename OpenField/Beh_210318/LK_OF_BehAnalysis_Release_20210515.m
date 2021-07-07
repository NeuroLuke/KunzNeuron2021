%==========================================================================
% Behavioral analysis.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clc; clear; close all;
opengl hardware;

% set rng
rng(1);

% paths
paths           = [];
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.save      = paths.beh;
paths.saveSubj  = strcat(paths.beh, 'SubjectPlots\');

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

% preallocate
numTrials   	= nan(length(subjects), 1); % number of trials
expDuration  	= nan(length(subjects), 1); % experiment duration

%% loop through all subjects
for iSub = 1:size(subjects, 1)
    
    %% data
    
    % trial information
    behFile   	= dir(strcat(paths.beh, subjects{iSub, 1}, '\trials.mat'));
    tmp       	= load(strcat(behFile.folder, '\', behFile.name));
    trials     	= tmp.trials;
    
    % information about number of trials and experiment duration (duration
    % between first ITI onset and last grab onset)
    numTrials(iSub, 1)      = size(trials, 1);
    expDuration(iSub, 1)    = (max(trials(:, 8)) - min(trials(:, 3))) / 60; % (minutes)    
end

%% report number of trials

fprintf('Number of trials.\n');
fprintf('Range: %d - %d trials.\n', min(numTrials), max(numTrials));
fprintf('Mean +/- SEM: %.3f +/- %.3f trials.\n', mean(numTrials), std(numTrials) ./ sqrt(size(numTrials, 1)));

%% report experiment duration

fprintf('Experiment duration.\n');
fprintf('Range: %.3f - %.3f minutes.\n', min(expDuration), max(expDuration));
fprintf('Mean +/- SEM: %.3f +/- %.3f minutes.\n', mean(expDuration), std(expDuration) ./ sqrt(size(expDuration, 1)));
