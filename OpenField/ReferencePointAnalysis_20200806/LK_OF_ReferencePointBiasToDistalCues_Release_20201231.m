%==========================================================================
% This script examines whether reference points are biased towards distal
% cues.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param       = [];
param.myRNG = 777;
rng(param.myRNG);

% paths
paths       = [];
paths.ELRPD = 'E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\';
paths.arena = 'E:\OpenField\Arena\';
paths.save  = 'E:\OpenField\ReferencePointAnalysis_20200806\Bias2DistalCues_20201231\';
mkdir(paths.save);

%% ELRPD cell information

% load ELRPD cell analysis
ELRPD     	= load(strcat(paths.ELRPD, 'results.mat'), ...
    'allRes', 'subjects', 'locDir', 'direc');

% ELRPD cells
bELRPDCell  = cell2mat({ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELRPDCell));

% reference points
allCOMxy    = cell2mat({ELRPD.allRes(bELRPDCell).locdir_COMxy}');

%% test whether the allocentric directions of reference points are aligned with the distal cues

% report
fprintf('\nExamine tuning of allocentric directions of reference points towards distal landmarks.\n');

% allocentric directions of the reference points
alloDirs   	= atan2(allCOMxy(:, 2), allCOMxy(:, 1)); % atan2(y, x)

% information about landmarks
la         	= load(strcat(paths.arena, 'LandmarkAngles_20200609.mat'));

% number of preferred directions per landmark
Ntotal    	= size(alloDirs, 1); % total number of directions (one per ELRPD cell)
NperLM      = nan(size(la.LandmarkAngles, 1), 1); % number of directions per landmark
expNperLM   = nan(size(la.LandmarkAngles, 1), 1); % expected number of directions per landmark
for iLM = 1:size(NperLM, 1)
    NperLM(iLM, 1)      = sum(rad2deg(alloDirs) >= la.LandmarkAngles{iLM, 2} & ...
        rad2deg(alloDirs) < la.LandmarkAngles{iLM, 3});
    expNperLM(iLM, 1)   = Ntotal * range(cell2mat(la.LandmarkAngles(iLM, 2:3))) / 360;
end

% number of preferred directions per unique landmark
uniqueLM    = unique(la.LandmarkAngles(:, 1));
NperULM     = nan(size(uniqueLM, 1), 1);
expNperULM  = nan(size(uniqueLM, 1), 1);
for iULM = 1:size(uniqueLM)
    NperULM(iULM, 1)    = sum(NperLM(strcmp(la.LandmarkAngles(:, 1), uniqueLM{iULM})));
    expNperULM(iULM, 1) = sum(expNperLM(strcmp(la.LandmarkAngles(:, 1), uniqueLM{iULM})));
end

% round expected values to nearest integer
expNperULM  = round(expNperULM);

% perform chi-squared test to examine whether the number of allocentric
% directions is increased towards one of the distal cues
contT               = [NperULM'; expNperULM'];
[~, pChiS, X2ChiS]  = chi2cont(contT);
fprintf('Deviation of observed landmark-related directions from expected values: X2 = %.3f, p = %.3f.\n', ...
    X2ChiS, pChiS);
