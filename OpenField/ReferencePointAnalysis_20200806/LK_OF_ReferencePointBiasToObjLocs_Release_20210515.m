%==========================================================================
% This script examines whether reference points are biased towards object
% locations.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath('E:\OpenField\Functions\');

% settings
param               = [];
param.myRNG         = 1;
param.objectType    = 'objLocs'; % objLocs | resLocs
param.sliceSize     = 450; % for circular arc procedure
param.numSurrogates = 1001; % number of surrogates
rng(param.myRNG);

% paths
paths       = [];
paths.beh   = 'E:\OpenField\Beh_210318\';
paths.save  = 'E:\OpenField\ReferencePointAnalysis_20200806\Bias2ObjLocs_20201231\';
mkdir(paths.save);

%% ELRPD cell results

% load ELRPD cell analysis
ELRPD       = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'subjects', 'locDir');

% subjects
subjects    = ELRPD.subjects;

% ELRPD cells
bELRPDCell  = cell2mat({ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELRPDCell));
allUnitIdx  = cell2mat({ELRPD.allRes.idx}');

% reference points and their distance towards the center
allCOMxy    = cell2mat({ELRPD.allRes.locdir_COMxy}');
allCOMD2Ctr = pdist2(allCOMxy, [0, 0]);

%% object/response locations for each subject and session

% select locations per subject (i.e., object/response locations)
selLocs = cell(size(subjects, 1), 1);
fprintf('\nUnder examination: %s.\n', param.objectType);

% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    thisSubject = subjects{iSub};
    fprintf('\tProcessing behavioral data for %s.\n', thisSubject);
    
    % load behavioral information
    tmp     = load(strcat(paths.beh, thisSubject, filesep, 'trials.mat'));
    trials  = tmp.trials; % trialidx, objectidx, ITI, cue, retrieval, feedback, reencoding, grab, correctXY, dropXY, drop error
    
    % correct object locations
    objLocs = sortrows(unique(trials(~isnan(trials(:, 9)), [2, 9:10]), 'rows', 'stable'));
    
    % average response locations
    resLocs     = nan(size(objLocs, 1), 3);
    for iObj = 1:size(objLocs, 1)
        
        % response locations
        thisResLocs         = trials(trials(:, 2) == iObj - 1, [2, 11:12]);
        resLocs(iObj, :)    = nanmean(thisResLocs);
    end    
    
    % choose object locations or response locations
    if strcmp(param.objectType, 'objLocs')
        selLocs{iSub, 1}    = objLocs;
    elseif strcmp(param.objectType, 'resLocs')
        selLocs{iSub, 1}    = resLocs;
    end
end

%% statistical assessment of the question whether reference points are significantly close to object locations using a circular-arc procedure

% report
fprintf('\nStatistical evaluation using a circular arc to create surrogate reference points.\n');

% create random locations within the arena to draw surrogate reference
% points
cfg         = [];
cfg.maxR    = 6000; % necessarily larger than the arena radius, e.g., 6000
cfg.minR    = 0;
cfg.N       = 2000001; % number of random locations
cfg.centerX = 0;
cfg.centerY = 0;
randLocs    = LK_RandomPointsInCircle_101119(cfg);

% distance of random locations to the center
allRandLocsD2Ctr    = pdist2(randLocs, [0, 0]);

% preallocate
minDEmp                 = nan(size(allCOMxy, 1), 1);
minDSurroAll            = cell(size(allCOMxy, 1), 1);
minDSurroSlice          = cell(size(allCOMxy, 1), 1);
minDSliceRanksSmaller  	= nan(size(allCOMxy, 1), 1); % ranks for empirical < surrogates (left-sided test)

% loop through cells
for iCell = 1:size(ELRPD.allRes, 1)
    
    % this cell's reference point
    thisCOMxy   = allCOMxy(iCell, :);
    
    % this subject's object locations
    thisLocs    = selLocs{allUnitIdx(iCell, 1)};
    
    % distance of reference point to object locations
    tmpD                = pdist2(thisLocs(:, 2:3), thisCOMxy);
    minDEmp(iCell, 1)   = min(tmpD);
        
    % use circular slice for comparison
    sliceBorders        = [allCOMD2Ctr(iCell, 1) - param.sliceSize, allCOMD2Ctr(iCell, 1) + param.sliceSize];
    logIdx              = allRandLocsD2Ctr >= sliceBorders(1) & allRandLocsD2Ctr < sliceBorders(2);
    
    % distances of surrogate reference points to object locations
    tmpD                        = pdist2(thisLocs(:, 2:3), randLocs);
    minDSurroAll{iCell, 1}      = min(tmpD, [], 1); % minimum across objects (first dimension)
    minDSurroSlice{iCell, 1}    = min(tmpD(:, logIdx), [], 1);
    
    % estimate how often the empirical distance is smaller than the
    % surrogate distances
    minDSliceRanksSmaller(iCell, 1) = sum(minDEmp(iCell, 1) < minDSurroSlice{iCell, 1}) / sum(~isnan(minDSurroSlice{iCell, 1}));
end

% test ranks against chance level to detect general bias
[~, p, ~, stats]    = ttest(minDSliceRanksSmaller(bELRPDCell), 0.5);
fprintf('Are reference points significantly close/far from object locations? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% identify ELRPD cells with significantly close reference points
fprintf('Number of ELRPD cells with a reference point significantly close to the closest object location: %d.\n', sum(bELRPDCell & minDSliceRanksSmaller > 0.95));

%% illustration of how the surrogate minimum distances are created - estimation

% example cell
exCellIdx   = 679; % index of example cell
exMinDSurro = transpose(minDSurroAll{exCellIdx, 1});

% distribute minimum distances onto grid
xEdges      = -5000:50:5000;
xCenters    = movmean(xEdges, 2, 'endpoints', 'discard');
yEdges      = -5000:50:5000;
yCenters    = movmean(yEdges, 2, 'endpoints', 'discard');
M           = nan(numel(yEdges), numel(xEdges)); % average minimum distance at each location
for iX = 1:numel(xEdges) - 1
    parfor iY = 1:numel(yEdges) - 1
        if sqrt(xCenters(iX) ^ 2 + yCenters(iY) ^ 2) > 5000
            continue;
        end
        logIdx      = randLocs(:, 1) >= xEdges(iX) & randLocs(:, 1) < xEdges(iX + 1) & randLocs(:, 2) >= yEdges(iY) & randLocs(:, 2) < yEdges(iY + 1);
        M(iY, iX)   = mean(exMinDSurro(logIdx));
    end
end

%% illustration of how the surrogate minimum distances are created - figure

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 9, 7]);
ax1 = axes('units', 'centimeters', 'position', [0.6, 0.6, 6.5, 6.5], 'visible', 'off');
hold on;
% minimum distances and object locations
myImag = imagesc(xCenters, yCenters, M, ...
    'AlphaData', ~isnan(M));
plot(selLocs{allUnitIdx(exCellIdx, 1)}(:, 2), selLocs{allUnitIdx(exCellIdx, 1)}(:, 3), '.', ...
    'Color', [0, 0, 0], 'MarkerSize', 12);
% reference point
plot(allCOMxy(exCellIdx, 1), allCOMxy(exCellIdx, 2), 'o', ...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 10, 'LineWidth', 2);
% boundary
tmpX = cos(0:0.001:2*pi);
tmpY = sin(0:0.001:2*pi);
patch([tmpX .* 5000, fliplr(tmpX) .* 6250], [tmpY .* 5000, fliplr(tmpY) .* 6250], [1, 1, 1], ...
    'EdgeColor', 'none');
plot(5000 .* tmpX, 5000 .* tmpY, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% slice for ranking
plot((allCOMD2Ctr(exCellIdx, 1) + param.sliceSize) .* cos(0:0.085:2*pi), (allCOMD2Ctr(exCellIdx, 1) + param.sliceSize) .* sin(0:0.085:2*pi), '.', ...
    'Color', [0.75, 0.75, 0.75], 'LineWidth', 2, 'MarkerSize', 8); % outer slice border
plot((allCOMD2Ctr(exCellIdx, 1) - param.sliceSize) .* cos(0:0.1:2*pi), (allCOMD2Ctr(exCellIdx, 1) - param.sliceSize) .* sin(0:0.1:2*pi), '.', ...
    'Color', [0.75, 0.75, 0.75], 'LineWidth', 2, 'MarkerSize', 8); % inner slice border
colormap jet;
cb = colorbar;
set(cb, ...
    'units', 'centimeters', 'position', [6.8, 5.1, 0.4, 1.5], ...
    'limits', [0, max(cb.Limits)], 'ytick', [0, max(cb.Limits)], 'ticklabels', {'Close', 'Far'});
set(gca, ...
    'xlim', [-5400, 5400], 'ylim', [-5400, 5400], ...
    'ydir', 'reverse'); % flip y-axis
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
print(f, strcat(paths.save, 'MinD2', param.objectType, '_illustration'), '-dtiff', '-r450');

%--- illustration of how the surrogate minimum distances are created - histogram

% adjust colormap and data to display
myColors        = jet;
sliceBorders    = [allCOMD2Ctr(exCellIdx, 1) - param.sliceSize, allCOMD2Ctr(exCellIdx, 1) + param.sliceSize];
logIdx          = allRandLocsD2Ctr >= sliceBorders(1) & allRandLocsD2Ctr < sliceBorders(2);
thisD           = minDSurroAll{exCellIdx, 1}(1, logIdx); % object distances for this circular arc
[N, edges]      = histcounts(thisD, linspace(round(min(cb.Limits)), round(max(cb.Limits)), size(myColors, 1) + 1)); % histogram
x               = movmean(edges, 2, 'endpoints', 'discard');

% create figure for histogram
f = figure('units', 'centimeters', 'position', [2, 2, 5, 4]);
axes('units', 'centimeters', 'position', [0.9, 1, 2.5, 2.2]);
hold on;
for iX = 1:numel(x)
    if N(iX) < 1
        continue;
    end
    bar(x(iX), N(iX), ...
        'FaceColor', myColors(iX, :), 'EdgeColor', 'none', ...
        'BarWidth', unique(diff(edges)));    
end
xline(minDEmp(exCellIdx, :), ...
    'LineWidth', 2, 'Color', [0, 0, 0]);
t = text(0.175, 0.875, num2str(minDSliceRanksSmaller(exCellIdx), '%.2f'), ...
    'units', 'normalized');
tmpAx = get(gca);
set(gca, ...
    'xlim', [0, myceil(max(tmpAx.XLim), -1)], 'xtick', [0, myceil(max(tmpAx.XLim), -1)], 'xticklabel', {'Close', 'Far'}, 'xdir', 'reverse', ...
    'ytick', [], ...
    'tickdir', 'out', 'ticklength', [0.04, 0.04]);
yl = ylabel('Probability');
set([gca, yl, t], ...
    'fontunits', 'centimeters', 'fontsize', 0.45);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'MinD2', param.objectType, '_EmpiricalVSSurrogates'), '-dtiff', '-r450');

%% histogram of ELRPD-cell ranks regarding the minimum distance of reference points towards the closest object

% settings
xSpacing    = 0.05;
xEdges      = 0:xSpacing:1;
xMeans      = movmean(xEdges, 2, 'endpoints', 'discard');

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [1, 2.2468, 5.6273, 5.1443]);
hold on;
myH = histogram(minDSliceRanksSmaller(bELRPDCell), xEdges, ...
    'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
xSig = xMeans(xMeans > 0.95);
ySig = myH.BinCounts(xMeans > 0.95);
bar(xSig, ySig, ...
    'BarWidth', 0.05, 'EdgeColor', 'none', ...
    'FaceColor', [1, 0, 0]);
% cutoff for significant cells
xline(0.95, '-', ...
    'LineWidth', 2, 'Color', [0, 0, 0]);
% chance level for each bin
yline(sum(bELRPDCell) / myH.NumBins, '--', ...
    'LineWidth', 1);
xl = xlabel({'Proximity to closest', 'object location (ranked)'});
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
tmpAx = get(gca);
set(gca, ...
    'xlim', [min(xEdges) - xSpacing, max(xEdges) + xSpacing], ...
    'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ELRPDCells_MinD2', param.objectType, 'Ranks_20210406'), '-dtiff', '-r450');
