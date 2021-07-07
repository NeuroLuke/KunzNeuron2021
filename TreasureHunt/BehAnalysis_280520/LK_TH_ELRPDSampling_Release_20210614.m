%==========================================================================
% Analysis of sampling of egocentric directions towards local reference
% points (LRPs).
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));
rng(111);

% paths
paths           = [];
paths.beh       = 'E:\TreasureHunt\Beh_210420\';
paths.save      = 'E:\TreasureHunt\BehAnalysis_280520\';

% settings
param           = [];
param.diameter  = 110;
param.XLim      = [315, 425];
param.ZLim      = [305, 415];

% load behavioral data from ELRPD-cell analysis
r   = load('E:\TreasureHunt\ELRPDCellAnalysis_220420\20200722_Dir30_Loc6x6_EgoLocDir30\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', ...
    'locDir', 'subjects', 'allRes', 'allBeh', 'param');

% unit indices
allUnitIdx   	= cell2mat({r.allRes.idx}');

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
for iLRP = 1:size(r.locDir.LRPs, 1)
    
    % skip LRPs that were excluded from the analysis
    if r.locDir.LRPs2Exclude(iLRP, 1) == true
        continue;
    end
    
    % for each subject, get distribution of egocentric direcions towards
    % this local reference point
    allProb = [];
    behIdx  = 1;
    for iSub = 1:size(r.subjects, 1)
        
        % session indices
        sessIdx = unique(allUnitIdx(allUnitIdx(:, 1) == iSub, 2));
        
        % loop through sessions
        for iSess = 1:size(sessIdx, 1)
                        
            % get dwell time probability for each angular bin
            thisSessIdx = sessIdx(iSess);
            firstIdx    = find(allUnitIdx(:, 1) == iSub & allUnitIdx(:, 2) == thisSessIdx, 1, 'first');
            
            % get data from timepoints that entered the analysis
            bMask4Analysis  = r.allRes(firstIdx).bMask4Analysis;
            thisProb        = histcounts(r.allBeh(behIdx).egoLocDir{iLRP}(bMask4Analysis), r.locDir.angleEdges, ...
                'normalization', 'probability');
            allProb         = cat(1, allProb, thisProb);
            
            % go to next behavioral session
            behIdx          = behIdx + 1;
        end
    end
    
    % mean and STD across subjects
    th  = r.locDir.angleCenters; % angles
    m   = mean(allProb); % mean
    ste = std(allProb); % std
    
    % average distribution of egocentric directions towards this LRP
    axes('units', 'normalized', 'position', [0.88 * (1 - (r.locDir.LRPs(iLRP, 1) - min(param.XLim)) / param.diameter) + 0.05, ...
        0.88 * (1 - (r.locDir.LRPs(iLRP, 2) - min(param.ZLim)) / param.diameter) + 0.045, 0.065, 0.0625]);
    hold on;
    patch([th, fliplr(th)], [m + ste, fliplr(m - ste)], [0.7, 0.7, 0.7], 'EdgeColor', 'none');
    plot(th, m, 'Color', [0, 0, 0], 'LineWidth', 1);
    hold off;
    set(gca, ...
        'xlim', [-pi, pi], 'ylim', [0.05, 0.13], ...
        'xdir', 'reverse', ... % because negative angles mean "to the right"
        'fontunits', 'centimeters', 'fontsize', 0.17, ...
        'ticklength', [0.05, 0.05]);
    % add xlabel
    if r.locDir.LRPs2Exclude(iLRP - 1) == true || r.locDir.LRPs(iLRP, 2) == max(r.locDir.LRPs(:, 2))
        set(gca, 'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {'', 'R', 'A', 'L', 'B'}); % because x-axis is reversed
    else
        set(gca, 'xtick', []);
    end
    % add ylabel
    if r.locDir.LRPs(iLRP, 1) == max(r.locDir.LRPs(:, 1)) || r.locDir.LRPs2Exclude(iLRP - r.locDir.locRes) == true
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
