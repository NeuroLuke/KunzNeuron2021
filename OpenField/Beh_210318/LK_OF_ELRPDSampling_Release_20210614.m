%==========================================================================
% Analysis of sampling of egocentric directions towards local reference
% points (LRPs).
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param         	= [];
param.radius   	= 4950;
param.diameter	= param.radius * 2;
param.myRNG   	= 111;
rng(param.myRNG);

% paths
paths       = [];
paths.beh   = 'E:\OpenField\Beh_210318\';
paths.save  = paths.beh;

% load behavioral data from ELRPD-cell analysis
r   = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'locDir', 'subjects', 'allRes', 'allBehLRP', 'param');

% identify unique sessions
allUnitIdx          = cell2mat({r.allRes.idx}');
[~, idxUniqueSess]  = unique(allUnitIdx(:, 1));

%% ELRPD sampling as a function of allocentric location

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 13.5, 13.5]);

% enhance axes
axes('units', 'normalized', 'position', [0, 0, 1, 1], ...
    'visible', 'off');
text(0.475, 0.015, 'Egocentric bearing', 'units', 'normalized', ...
    'fontunits', 'centimeters', 'fontsize', 0.25);
text(0.01, 0.475, 'Probability', 'units', 'normalized', ...
    'Rotation', 90, ...
    'fontunits', 'centimeters', 'fontsize', 0.25);

% loop trough local reference points
th  = cell(size(r.locDir.LRPs, 1), 1); % angles
m   = cell(size(r.locDir.LRPs, 1), 1); % mean
ste = cell(size(r.locDir.LRPs, 1), 1); % STD
for iLRP = 1:size(r.locDir.LRPs, 1)
    
    % skip LRPs that were excluded from the analysis
    if r.locDir.LRPs2Exclude(iLRP, 1) == true
        continue;
    end
    
    % for each session, distribution of egocentric directions towards this
    % local reference point
    prob    = cell(size(r.subjects, 1), 1);
    for iSub = 1:size(r.subjects, 1)
        
        % get dwell time probability for each angular bin
        bMask4Analysis  = r.allRes(idxUniqueSess(iSub, 1)).bMask4Analysis;
        prob{iSub, 1}   = histcounts(r.allBehLRP{iSub, iLRP}.egoLocDir(bMask4Analysis), r.locDir.angleEdges, ...
            'normalization', 'probability');
    end
    
    % egocentric angles, mean, and STD
    th{iLRP}    = r.locDir.angleCenters';
    m{iLRP}  	= mean(cell2mat(prob));
    ste{iLRP}   = std(cell2mat(prob));
        
    % average distribution of egocentric direction towards this LRP
    axes('units', 'normalized', 'position', [0.88 * (r.locDir.LRPs(iLRP, 1) + param.radius) / param.diameter + 0.05, ...
        0.88 * (1 - (r.locDir.LRPs(iLRP, 2) + param.radius) / param.diameter) + 0.045, 0.065, 0.0625]);
    hold on;
    patch([th{iLRP}, fliplr(th{iLRP})], [m{iLRP} + ste{iLRP}, fliplr(m{iLRP} - ste{iLRP})], [0.7, 0.7, 0.7], ...
        'EdgeColor', 'none');
    plot(th{iLRP}, m{iLRP}, ...
        'Color', [0, 0, 0], 'LineWidth', 1);
    hold off;
    set(gca, ...
        'xlim', [-pi, pi], ...
        'ylim', [0.05, 0.13], ...
        'fontunits', 'centimeters', 'fontsize', 0.17, ...
        'ticklength', [0.05, 0.05]);
    % add xlabel
    if r.locDir.LRPs2Exclude(iLRP + 1) == true || r.locDir.LRPs(iLRP, 2) == max(r.locDir.LRPs(:, 2))
        set(gca, 'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {'B', 'L', 'A', 'R', ''});
    else
        set(gca, 'xtick', []);
    end
    % add ylabel
    if r.locDir.LRPs(iLRP, 1) == min(r.locDir.LRPs(:, 1)) || r.locDir.LRPs2Exclude(iLRP - 12) == true
        set(gca, 'ytick', [0.06, 0.12]);
    else
        set(gca, 'ytick', []);
    end
    % add title
    title([num2str(r.locDir.LRPs(iLRP, 1)), '/', num2str(r.locDir.LRPs(iLRP, 2))], ...
        'fontunits', 'centimeters', 'fontsize', 0.17, 'fontweight', 'normal', ...
        'units', 'normalized', 'position', [0.55, 0.85, 0]);
    drawnow;
end

% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'DistributionOfELRPDsDependingOnLRP_20210409'), '-dtiff', '-r450');
