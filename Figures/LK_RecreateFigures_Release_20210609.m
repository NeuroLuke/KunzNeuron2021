%==========================================================================
% This script recreates the figures.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('C:\GitCode\KunzNeuron2021\OpenField\Functions\'));
addpath(genpath('C:\GitCode\KunzNeuron2021\TreasureHunt\Functions\'));

% paths
paths       = [];
paths.data  = 'C:\GitCode\KunzNeuron2021\Figures\';

%% Figure 1C

% data
dt  = load(strcat(paths.data, 'Fig_1C_data.mat'));

% histogram of all memory accuracies
figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
hold on;
myH = histogram(cell2mat(dt.allMemAcc), 0:0.05:1, ...
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
drawnow;

%% Figure 1D

% data
dt  = load(strcat(paths.data, 'Fig_1D_data.mat'));

% memory performance in the last vs. first trial
figure('units', 'centimeters', 'position', [5, 5, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.55, 3.8, 4]);
hold on;
for iSub = 1:size(dt.allMemAcc_trialFirstLast, 1)
    plot([1, 2], dt.allMemAcc_trialFirstLast(iSub, :), '-', ...
        'Color', rgb('gray'));
end
plot([1, 2], mean(dt.allMemAcc_trialFirstLast), '-', ...
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
drawnow;

%% Figure 2, A-E

% loop through the different subpanels
panels  = {'A'; 'B'; 'C'; 'D'; 'E'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt          = load(strcat(paths.data, 'Fig_2', panels{iPanel}, '_data.mat'));
    dt          = dt.dt;
    dt.visible  = 'on';
    
    % ELRPD tuning
    LK_OF_PlotELRPDModulation_20210408(dt);
    drawnow;
end

%% Figure 3A

% data
dt  = load(strcat(paths.data, 'Fig_3A_data.mat'));

% create Venn diagram
figure('units', 'centimeters', 'position', [2, 5, 12, 8]);
venn(dt.A, dt.I, 'ErrMinMode', 'None', ...
    'FaceColor', {[0, 0.5, 1], [0, 0, 1], [0, 1, 0]}, ...
    'LineWidth', 2);
axis equal tight off;
drawnow;

%% Figure 3B

% data
dt  = load(strcat(paths.data, 'Fig_3B_data.mat'));

% create bar plot for egocentric bearing cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Egocentric bearing cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 3C

% data
dt  = load(strcat(paths.data, 'Fig_3C_data.mat'));

% plot reference points
figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
hold on;
plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 5);
plot(cos(0:0.001:(2*pi)) .* dt.borderCtrPeri, sin(0:0.001:(2*pi)) .* dt.borderCtrPeri, ':', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
scatter(dt.ELRPDCells_allCOMxy(dt.ELRPDCells_allCOMD2Ctr_Group == 1, 1), ...
    dt.ELRPDCells_allCOMxy(dt.ELRPDCells_allCOMD2Ctr_Group == 1, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % center reference points
scatter(dt.ELRPDCells_allCOMxy(dt.ELRPDCells_allCOMD2Ctr_Group == 2, 1), ...
    dt.ELRPDCells_allCOMxy(dt.ELRPDCells_allCOMD2Ctr_Group == 2, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % periphery reference points
hold off;
set(gca, ...
    'ydir', 'reverse', ... % negative y values at top
    'xlim', [-5400, 5400], 'xtick', [-4950, 4950], ...
    'ylim', [-5400, 5400], 'ytick', [-4950, 4950], ...
    'tickdir', 'out', 'ticklength', [0, 0], 'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
ytickangle(90);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 3D

% data
dt  = load(strcat(paths.data, 'Fig_3D_data'));

% histogram for distance of reference points (COMs) to the environmental
% center
figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.55, 1.75, 6, 6]);
hold on;
% histograms of both groups
H1 = histogram(dt.ELRPDCells_allCOMD2Ctr(dt.ELRPDCells_allCOMD2Ctr_Group == 1), dt.rEdges, ...
    'facecolor', rgb('darkgreen'), 'edgecolor', 'none', 'facealpha', 1); % center reference points
H2 = histogram(dt.ELRPDCells_allCOMD2Ctr(dt.ELRPDCells_allCOMD2Ctr_Group == 2), dt.rEdges, ...
    'facecolor', rgb('limegreen'), 'edgecolor', 'none', 'facealpha', 1); % periphery reference points
H3 = stairs(dt.rEdges, [dt.numCOMsPerRExp; dt.numCOMsPerRExp(end)], ...
    'Color', [0.5, 0.5, 0.5]); % expected distribution
axis square;
set(gca, ...
    'xlim', [min(dt.rEdges) - 100, max(dt.rEdges) + 100], 'xtick', [0, 2500, 5000], ...
    'ylim', [0, max([H1.BinCounts, H2.BinCounts, H3.YData]) + 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
hold off;
xl = xlabel('Distance to center (vu)');
yl = ylabel('Count');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5)
uistack(gca, 'top');
drawnow;

%% Figure 3E

% loop through the different subpanels
panels  = {'Left'; 'Right'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_3E', panels{iPanel}, '_data.mat'));
    
    % create polarhistogram
    figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
    axis off;
    
    % angle-circle
    ax1 = axes('units', 'centimeters', 'position', [1.5, 1.5, 5.5, 5.5]);
    hold on;
    % plot angle lines
    for iAng = 0:30:165
        plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
            'Color', rgb('gray'), 'LineWidth', 1);
    end
    % plot circle
    tmpx = (-pi:0.001:pi) + dt.rotAngle; % rotate so that red is pointing up
    tmpc = linspace(0, 1, length(tmpx));
    patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
        [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
        'EdgeColor', 'interp');
    % plot angle text
    t1 = text(1.125, 0, 'L'); % L here, because x-axis will be reversed
    t2 = text(0, 1.15, 'A');
    t3 = text(-1.125, 0, 'R'); % R here, because x-axis will be reversed
    t4 = text(0, -1.15, 'B');
    set([t1, t2, t3, t4], ...
        'FontUnits', 'centimeters', 'FontSize', dt.myFontSize, ...
        'HorizontalAlignment', 'center');
    % indicate maximum
    text(0.175, 0.15, ['max ', num2str(max(dt.NcircM), '%.0f')], ...
        'units', 'normalized', ...
        'horizontalalignment', 'right', 'verticalalignment', 'top', ...
        'fontunits', 'centimeters', 'fontsize', dt.myFontSize)
    hold off;
    % reverse x-axis so that -90° (i.e., left) points to the left
    set(gca, ...
        'xdir', 'reverse', ...
        'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
    myCM = circshift(hsv, size(hsv, 1) / 2); % so that 0° is red
    colormap(ax1, myCM);
    axis off;
    
    % observations per angle
    axes('units', 'centimeters', 'position', ax1.Position);
    maxN = max(dt.NcircM);
    hold on;
    for iA = 1:size(dt.angleCenters, 1)
        tmpAngles   = transpose([0, dt.angleCenters(iA) - deg2rad(dt.angularRes / 2):...
            0.001:dt.angleCenters(iA) + deg2rad(dt.angularRes / 2), 0] + dt.rotAngle);
        tmpRadii    = [0; repmat(dt.NcircM(iA), numel(tmpAngles) - 2, 1); 0] ./ maxN;
        [x, y]      = pol2cart(tmpAngles, tmpRadii);
        patch(x, y, [0.7, 0.7, 0.7], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(x, y, '-', ...
            'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    end
    hold off;
    % adjust axes
    set(gca, ...
        'xdir', 'reverse', ...
        'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
    axis off;
    drawnow;
end

%% Figure 3F

% loop through the different subpanels
panels  = {'Left'; 'Right'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt          = load(strcat(paths.data, 'Fig_3F', panels{iPanel}, '_data.mat'));
    dt          = dt.ld;
    dt.visible  = 'on';
    
    % create figure for linear distance tuning
    LK_PlotLinearDistModulation_20201011(dt);
    drawnow;
end

%% Figure 3G

% data
dt  = load(strcat(paths.data, 'Fig_3G_data.mat'));

% plot all tuning curves
figure('units', 'centimeters', 'position', [2, 2, 7, 16]);
axes('units', 'centimeters', 'position', [1, 1.75, 5, 12.75]);
imagesc(dt.binCenters, 1:size(dt.allFR, 1), dt.allFR(dt.sortIdx, :), ...
    'alphadata', dt.alphaData(dt.sortIdx, :));
cb = colorbar('southoutside');
cb.Units = 'centimeters';
cb.Position = [2.25, 14.75, 2.5, 0.3];
ylabel(cb, 'FR', ...
    'units', 'normalized', 'position', [0.5, 1.55]);
cb.Ticks = [0, 1];
cb.TickLabels = {'min', 'max'};
xl = xlabel('Distance (vu)');
yl = ylabel('Cell index', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set(gca, ...
    'xtick', [min(dt.binCenters), max(dt.binCenters)], ...
    'ytick', [min(dt.sortIdx), max(dt.sortIdx)], ...
    'ticklength', [0, 0]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
colormap jet;
drawnow;

%% Figure 3, H-J

% loop through different subpanels
panels  = {'H'; 'I'; 'JLeft'; 'JRight'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt          = load(strcat(paths.data, 'Fig_3', panels{iPanel}, '_data.mat'));
    dt          = dt.bd;
    dt.visible  = 'on';
    
    % plot bearing-distance map of firing rates
    LK_PlotBearingDistanceMap_20201011(dt);
    drawnow;
end

%% Figure 3K

% data
dt  = load(strcat(paths.data, 'Fig_3K_data.mat'));

% create figure for bearing-distance extent of bearing-distance fields
figure('units', 'centimeters', 'position', [2, 2, 7, 7]);
hold on;
plot(100 .* dt.xyExtentFieldPerc(:, 1), 100 .* dt.xyExtentFieldPerc(:, 2), 'o', ...
    'Color', [0.7, 0.7, 0.7]); % in percent
plot(100 .* dt.xyExtentFieldPerc(dt.bDistCellField, 1), 100 .* dt.xyExtentFieldPerc(dt.bDistCellField, 2), 'o', ...
    'Color', [0.7, 0.7, 0.7], 'MarkerFaceColor', rgb('green')); % in percent
hold off;
xl = xlabel('Bearing extent (%)');
yl = ylabel('Distance extent (%)');
axis equal;
set(gca, ...
    'xlim', [0, 100], 'xtick', linspace(0, 100, 5), 'ylim', [0, 100], 'ytick', linspace(0, 100, 5), ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);

% create figure for total size of bearing-distance fields
figure('units', 'centimeters', 'position', [2, 2, 4, 4]);
axes('units', 'centimeters', 'position', [1, 1.6, 2.5, 2.1]);
hold on;
histogram(100 .* dt.extentFieldPerc(~dt.bDistCellField), dt.xEdges, ...
    'FaceColor', [0.7, 0.7, 0.7], 'facealpha', 0.5);
histogram(100 .* dt.extentFieldPerc(dt.bDistCellField), dt.xEdges, ...
    'FaceColor', rgb('green'), 'facealpha', 0.5);
xl = xlabel('Extent (%)');
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
tmpAx = get(gca);
set(gca, ...
    'xlim', [min(dt.xEdges), max(dt.xEdges)], 'xtick', [min(dt.xEdges), max(dt.xEdges)], ...
    'ytick', tmpAx.YLim, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 3L

% data
dt  = load(strcat(paths.data, 'Fig_3L_data.mat'));

% summary 2D bearing-distance plot
figure('units', 'centimeters', 'position', [2, 2, 7.5, 8], 'Color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [1.5, 1.75, 4, 4]);
hold on;
imagesc(dt.bearDist.bearingBinCenters, dt.bearDist.distanceBinCenters, dt.sumBearingDistanceField, ...
    'alphadata', dt.sumBearingDistanceField ~= 0);
hold off;
myCM = jet;
colormap(myCM);
cb = colorbar;
ylabel(cb, 'Count', 'Rotation', 0, ...
    'units', 'normalized', 'position', [1.75, 0.5], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
cb.Units = 'centimeters';
cb.Position = [5.8, 6.05, 0.25, 1.25];
cb.Limits = round([1, max(cb.Limits)], 1);
cb.Ticks = round([1, max(cb.Limits)], 1);
cb.Label.FontUnits = 'centimeters';
cb.Label.FontSize = 0.4;
xl = xlabel('Bearing');
yl = ylabel('Distance (vu)', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
idxValidY = find(sum(dt.sumBearingDistanceField, 2) > 0, 1, 'first'):find(sum(dt.sumBearingDistanceField, 2) > 0, 1, 'last') + 1;
set(gca, ...
    'xlim', [min(dt.bearDist.bearingBinEdges), max(dt.bearDist.bearingBinEdges)], 'xtick', linspace(min(dt.bearDist.bearingBinEdges), max(dt.bearDist.bearingBinEdges), 5), ...
    'xticklabel', {'B', 'L', 'A', 'R', ''}, ...
    'ylim', [0, max(dt.bearDist.distanceBinEdges(idxValidY))], 'ytick', [min(dt.bearDist.distanceBinEdges), max(dt.bearDist.distanceBinEdges(idxValidY))], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'Fontunits', 'centimeters', 'Fontsize', 0.5);

% distance summary
ax2 = axes('units', 'centimeters', 'position', [5.8, 1.75, 1.25, 4]);
sumDist = sum(dt.sumBearingDistanceField, 2);
plot(sumDist(idxValidY), dt.bearDist.distanceBinCenters(idxValidY), '-', ...
    'Color', [0.5, 0.5, 0.5]);
set(gca, ...
    'xlim', [0, max(sumDist(idxValidY))], 'xtick', [], ...
    'ylim', [0, max(dt.bearDist.distanceBinEdges(idxValidY))], 'ytick', [], ...
    'box', 'off');

% bearing summary
ax3 = axes('units', 'centimeters', 'position', [1.5, 6.05, 4, 1.25]);
sumBearing = sum(dt.sumBearingDistanceField, 1);
plot(dt.bearDist.bearingBinCenters, sumBearing, '-', ...
    'Color', [0.5, 0.5, 0.5]);
set(gca, ...
    'xlim', [min(dt.bearDist.bearingBinEdges), max(dt.bearDist.bearingBinEdges)], 'xtick', [], ...
    'ylim', [0, max(sumBearing)], 'ytick', [], ...
    'ytick', [], ...
    'box', 'off');
drawnow;

%% Figure 3M

% data
dt  = load(strcat(paths.data, 'Fig_3M_data.mat'));

% create bar plot for memory cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Memory cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 3N

% loop through different subpanels
panels  = {'Left'; 'Right'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt          = load(strcat(paths.data, 'Fig_3N', panels{iPanel}, '_data.mat'));
    dt          = dt.dt;
    dt.visible  = 'on';
    
    % plot relationship between firing rate and memory performance
    LK_PlotMemoryTuning_20210515(dt);
    drawnow;
end

%% Figure 3O

% data
dt  = load(strcat(paths.data, 'Fig_3O_data.mat'));

% create bar plot
figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [2.9, 2.3, 3.6, 5.6]);
hold on;
barh(1:numel(dt.nameOfCellType), repmat(100, numel(dt.nameOfCellType), 1), 'FaceColor', [1, 1, 1]);
barh(1:numel(dt.nameOfCellType), dt.percPerCellType * 100, 'FaceColor', [0, 0, 0]);
set(gca, ...
    'xlim', [0, 100], ...
    'ylim', [0.4, size(dt.nameOfCellType, 1) + 0.6], ...
    'ytick', 1:size(dt.nameOfCellType, 1), 'yticklabel', '', ...
    'ydir', 'reverse', ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
myAxes = get(gca);
for iLabel = 1:size(dt.nameOfCellType, 1)
    text(-8, myAxes.YTick(iLabel), dt.nameOfCellType{iLabel}, 'HorizontalAlignment', 'right', 'fontunits', 'centimeters', 'fontsize', 0.5);
end
xl = xlabel({'% memory cells', 'for each cell type'});
set([gca, xl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 4A

% data
dt          = load(strcat(paths.data, 'Fig_4A_data.mat'));
dt          = dt.dt;
dt.visible  = 'on';

% plot object cell
LK_PlotFRObj_20210515(dt);
drawnow;

%% Figure 4B

% data
dt  = load(strcat(paths.data, 'Fig_4B_data.mat'));

% distance between preferred objects
figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [1, 2.24682910989324, 5.25, 5.14430630677342]);
hold on;
myH = histogram(dt.meanD_prefObjLocsSurro, 2500:100:4500, ...
    'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
plot([dt.meanD_prefObjLocs, dt.meanD_prefObjLocs], [0, max(myH.Values)], ...
    'Color', [1, 0, 0], 'LineWidth', 2);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [2500, 4500], 'xtick', [2500, 4500], ...
    'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel({'Distance between', 'pref. objects (vu)'});
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);

% number of preferred objects
figure('units', 'centimeters', 'position', [2, 2, 4, 4]);
h1 = histogram(dt.objCells_numPrefObj, ...
    'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', [1, 1, 1]);
tmpAx = get(gca);
set(gca, ...
    'xtick', movmean(h1.BinEdges, 2, 'endpoints', 'discard'), ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('# pref. objects');
yl = ylabel('Cell count');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.4);
drawnow;

%% Figure 4C

% data
dt  = load(strcat(paths.data, 'Fig_4C_data.mat'));

% create bar plot for object cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Object cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 4D

% data
dt  = load(strcat(paths.data, 'Fig_4D_data.mat'));

% create Venn diagram
figure('units', 'centimeters', 'position', [2, 5, 12, 8]);
venn(dt.A, dt.I, ...
    'FaceColor', {rgb('orange'), [0.5, 1, 0]}, 'LineWidth', 2);
axis equal tight off;
drawnow;

%% Figure 4E

% data
dt  = load(strcat(paths.data, 'Fig_4E_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 9, 7], 'color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [0.6, 0.6, 6.5, 6.5], 'visible', 'off');
hold on;
% minimum distances and object locations
myImag = imagesc(dt.xCenters, dt.yCenters, dt.M, ...
    'AlphaData', ~isnan(dt.M));
plot(dt.exCell_selLocs(:, 2), dt.exCell_selLocs(:, 3), '.', ...
    'Color', [0, 0, 0], 'MarkerSize', 12);
% reference point
plot(dt.exCell_COMxy(1), dt.exCell_COMxy(2), 'o', ...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 10, 'LineWidth', 2);
% boundary
tmpX = cos(0:0.001:2*pi);
tmpY = sin(0:0.001:2*pi);
plot(5000 .* tmpX, 5000 .* tmpY, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% slice for ranking
plot((dt.exCell_COMD2Ctr + dt.param.sliceSize) .* cos(0:0.085:2*pi), (dt.exCell_COMD2Ctr + dt.param.sliceSize) .* sin(0:0.085:2*pi), '.', ...
    'Color', [0.75, 0.75, 0.75], 'LineWidth', 2, 'MarkerSize', 8); % outer slice border
plot((dt.exCell_COMD2Ctr - dt.param.sliceSize) .* cos(0:0.1:2*pi), (dt.exCell_COMD2Ctr - dt.param.sliceSize) .* sin(0:0.1:2*pi), '.', ...
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

% adjust colormap and data to display
myColors 	= jet;
thisD     	= dt.exCell_minDSurro; % object distances for this circular arc
[N, edges]  = histcounts(thisD, linspace(round(min(cb.Limits)), round(max(cb.Limits)), size(myColors, 1) + 1)); % histogram
x         	= movmean(edges, 2, 'endpoints', 'discard');

% create figure for histogram
figure('units', 'centimeters', 'position', [2, 2, 5, 4]);
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
xline(dt.exCell_minDEmp, ...
    'LineWidth', 2, 'Color', [0, 0, 0]);
tmpAx = get(gca);
set(gca, ...
    'xlim', [0, myceil(max(tmpAx.XLim), -1)], 'xtick', [0, myceil(max(tmpAx.XLim), -1)], 'xticklabel', {'Close', 'Far'}, 'xdir', 'reverse', ...
    'ytick', [], ...
    'tickdir', 'out', 'ticklength', [0.04, 0.04]);
yl = ylabel('Probability');
set([gca, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.45);
drawnow;

%% Figure 4F

% data
dt  = load(strcat(paths.data, 'Fig_4F_data.mat'));

% histogram of ranks
figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [1, 2.2468, 5.6273, 5.1443]);
hold on;
myH = histogram(dt.ELRPDCells_minDSliceRanksSmaller, dt.xEdges, ...
    'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
xSig = dt.xMeans(dt.xMeans > 0.95);
ySig = myH.BinCounts(dt.xMeans > 0.95);
bar(xSig, ySig, ...
    'BarWidth', 0.05, 'EdgeColor', 'none', ...
    'FaceColor', [1, 0, 0]);
% cutoff for significant cells
xline(0.95, '-', ...
    'LineWidth', 2, 'Color', [0, 0, 0]);
xl = xlabel({'Proximity to closest', 'object location (ranked)'});
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, xl, yl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
tmpAx = get(gca);
set(gca, ...
    'xlim', [min(dt.xEdges) - dt.xSpacing, max(dt.xEdges) + dt.xSpacing], ...
    'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
drawnow;

%% Figure 5B

% data
dt  = load(strcat(paths.data, 'Fig_5B_data.mat'));

% create figure for object-cell activity
myYLim  = [-1, 3.5];
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
patch([dt.time, fliplr(dt.time)], ...
    [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], rgb('orange'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3); % preferred
patch([dt.time, fliplr(dt.time)], ...
    [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], [0.7, 0.7, 0.7], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3); % unpreferred
plot(dt.time, dt.m{1}, '-', ...
    'Color', rgb('orange'), 'LineWidth', 2); % preferred
plot(dt.time, dt.m{2}, '-', ...
    'Color', [0.7, 0.7, 0.7], 'LineWidth', 2); % unpreferred
set(gca, ...
    'YLim', myYLim);
% add x- and y-axis
myxline(0, '--', [0, 0, 0], 'bottom');
myyline(0, '--', [0, 0, 0], 'bottom');
hold off;
set(gca, ...
    'xlim', [min(dt.time), max(dt.time)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Time (s)');
yl = ylabel('Firing rate (Hz)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 5, C-D

% loop through different panels
panels  = {'C'; 'D'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_5', panels{iPanel}, '_data.mat'));
    
    % create figure for preferred and unpreferred objects
    myYLim  = [-1, 3.5];
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], dt.myColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5); % object by egocentric bearing cells
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], dt.myColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % object by non-egocentric-bearing cells
    plot(dt.time, dt.m{1}, '-', ...
        'Color', dt.myColor, 'LineWidth', 2); % object by egocentric bearing cells
    plot(dt.time, dt.m{2}, ':', ...
        'Color', dt.myColor, 'LineWidth', 2); % object by non-egocentric-bearing cells
    set(gca, ...
        'YLim', myYLim);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(dt.time), max(dt.time)], ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure 5E

% data
dt  = load(strcat(paths.data, 'Fig_5E_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 8, 8], 'color', [1, 1, 1]);
ax1 = axes('units', 'centimeters', 'position', [0.6, 0.6, 6.5, 6.5], 'visible', 'off');
hold on;
% object locations and their distances to the reference point
for iObj = 1:size(dt.thisObjLocs, 1)
    if dt.bClose(iObj) == true
        myColor	= rgb('darkgreen'); % close objects
    else
        myColor = rgb('lightgreen'); % far objects
    end
    plot([dt.thisObjLocs(iObj, 1), dt.thisRP(1)], [dt.thisObjLocs(iObj, 2), dt.thisRP(2)], ':', ...
        'Color', myColor, 'LineWidth', 2);
end
plot(dt.thisObjLocs(~dt.bClose, 1), dt.thisObjLocs(~dt.bClose, 2), 'o', ...
    'MarkerSize', 10, 'Color', rgb('lightgreen'), 'MarkerFaceColor', rgb('lightgreen')); % far objects
plot(dt.thisObjLocs(dt.bClose, 1), dt.thisObjLocs(dt.bClose, 2), 'o', ...
    'MarkerSize', 10, 'Color', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen')); % close objects
% reference point
plot(dt.thisRP(1), dt.thisRP(2), 'o', ...
    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [1, 1, 1], 'MarkerSize', 10, 'LineWidth', 2);
% boundary
tmpX = cos(0:0.001:2*pi);
tmpY = sin(0:0.001:2*pi);
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
drawnow;

%% Figure 5F

% data
dt  = load(strcat(paths.data, 'Fig_5F_data.mat'));

% create figure
myYLim  = [-0.6, 0.8];
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
patch([dt.time, fliplr(dt.time)], ...
    [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], rgb('darkgreen'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
patch([dt.time, fliplr(dt.time)], ...
    [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], rgb('lightgreen'), ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(dt.time, dt.m{1}, '-', ...
    'Color', rgb('darkgreen'), 'LineWidth', 2);
plot(dt.time, dt.m{2}, '-', ...
    'Color', rgb('lightgreen'), 'LineWidth', 2);
set(gca, ...
    'YLim', myYLim, 'YTick', min(myYLim):0.2:max(myYLim));
% add x- and y-axis
myxline(0, '--', [0, 0, 0], 'bottom');
myyline(0, '--', [0, 0, 0], 'bottom');
hold off;
set(gca, ...
    'xlim', [min(dt.time), max(dt.time)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Time (s)');
yl = ylabel('Firing rate (Hz)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 6B, top and bottom

%- top: performance during (location-cued) object recall

% data
dt  = load(strcat(paths.data, 'Fig_6BTop_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 6, 5]);
axes('units', 'centimeters', 'position', [1.5, 2.2, 3.8, 2.5]);
hold on;
histogram(cell2mat(dt.objRecallPerf), 0:0.05:1, ...
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
drawnow;

%- bottom: performance during (object-cued) location recall

% data
dt  = load(strcat(paths.data, 'Fig_6BBottom_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 6, 5]);
axes('units', 'centimeters', 'position', [1.5, 2.2, 3.8, 2.5]);
hold on;
myH = histogram(dt.trialLocRecallPerf, 0:0.05:1, ...
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
drawnow;

%% Figure 6, C-G

% loop through the different subpanels
panels  = {'C'; 'D'; 'E'; 'F'; 'G'};
for iPanel = 1:size(panels, 1)
    
    % data for this panel
    dt          = load(strcat(paths.data, 'Fig_6', panels{iPanel}, '_data.mat'));
    dt          = dt.dt;
    dt.visible  = 'on';
    
    % ELRPD tuning
    LK_TH_PlotELRPDModulation_20210515(dt);
    drawnow;
end

%% Figure 6H

% data
dt  = load(strcat(paths.data, 'Fig_6H_data.mat'));

% create bar plot for egocentric bearing cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Egocentric bearing cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 6I

% data
dt  = load(strcat(paths.data, 'Fig_6I_data.mat'));

% plot reference points
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
hold on;
% boundary
plot(dt.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* dt.param.arenaRadius, dt.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* dt.param.arenaRadius, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 5);
% center/periphery boundary
plot(dt.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* dt.borderCtrPeri, dt.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* dt.borderCtrPeri, ':', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% reference points in center and periphery
scatter(dt.ELRPDCells_allCOMxz(dt.ELRPDCells_allCOMD2Ctr_Group == 1, 1), ...
    dt.ELRPDCells_allCOMxz(dt.ELRPDCells_allCOMD2Ctr_Group == 1, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1);
scatter(dt.ELRPDCells_allCOMxz(dt.ELRPDCells_allCOMD2Ctr_Group == 2, 1), ...
    dt.ELRPDCells_allCOMxz(dt.ELRPDCells_allCOMD2Ctr_Group == 2, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1);
hold off;
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ... % flip x- and y-axis
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], 'xtick', [min(dt.locDir.xCenters), max(dt.locDir.xCenters)], ...
    'ylim', [min(dt.locDir.zEdges), max(dt.locDir.zEdges)], 'ytick', [min(dt.locDir.zCenters), max(dt.locDir.zCenters)], ...
    'tickdir', 'out', 'ticklength', [0, 0], ...
    'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
ytickangle(90);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 6J

% data
dt  = load(strcat(paths.data, 'Fig_6J_data.mat'));

% histogram for distances of reference points (COMs) to the environmental
% center
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
axes('units', 'centimeters', 'position', [1.55, 1.75, 6, 6]);
hold on;
% histograms of both groups
H1 = histogram(dt.ELRPDCells_allCOMD2Ctr(dt.ELRPDCells_allCOMD2Ctr_Group == 1), dt.rEdges, ...
    'facecolor', rgb('darkgreen'), 'edgecolor', 'none', 'facealpha', 1); % center reference points
H2 = histogram(dt.ELRPDCells_allCOMD2Ctr(dt.ELRPDCells_allCOMD2Ctr_Group == 2), dt.rEdges, ...
    'facecolor', rgb('limegreen'), 'edgecolor', 'none', 'facealpha', 1); % periphery reference points
H3 = stairs(dt.rEdges, [dt.numCOMsPerRExp; dt.numCOMsPerRExp(end)], ...
    'Color', [0.5, 0.5, 0.5]);
axis square;
set(gca, ...
    'xlim', [min(dt.rEdges) - 1, max(dt.rEdges) + 1], 'xtick', [0, 25, 50], ...
    'ylim', [0, max([H1.BinCounts, H2.BinCounts, H3.YData]) + 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
hold off;
xl = xlabel('Distance to center (vu)');
yl = ylabel('Count');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure 6K

% loop through different panels
panels  = {'Top'; 'Bottom'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_6K', panels{iPanel}, '_data.mat'));
    
    % figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    axis off;
    % angle-circle
    ax1 = axes('units', 'centimeters', 'position', [1.5, 1.5, 5.5, 5.5]);
    hold on;
    % plot angle lines
    for iAng = 0:30:165
        plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
            'Color', rgb('gray'), 'LineWidth', 1);
    end
    % plot circle
    tmpx = (-pi:0.001:pi) + dt.rotAngle; % rotate so that red is pointing "north"
    tmpc = linspace(0, 1, length(tmpx));
    patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
        [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
        'EdgeColor', 'interp');
    % plot angle text
    t1 = text(1.125, 0, 'R');
    t2 = text(0, 1.15, 'A');
    t3 = text(-1.125, 0, 'L');
    t4 = text(0, -1.15, 'B');
    set([t1, t2, t3, t4], ...
        'FontUnits', 'centimeters', 'FontSize', dt.myFontSize, ...
        'HorizontalAlignment', 'center');
    % indicate maximum
    text(0.175, 0.15, ['max ', num2str(max(dt.NcircM), '%.0f')], ...
        'units', 'normalized', ...
        'horizontalalignment', 'right', 'verticalalignment', 'top', ...
        'fontunits', 'centimeters', 'fontsize', dt.myFontSize);
    hold off;
    % adjust x-axis
    set(gca, ...
        'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
    myCM    = flipud(circshift(hsv, size(hsv, 1) / 2)); % so that "ahead" is red and "right" is green
    colormap(ax1, myCM);
    axis off;
    
    % observations per angle
    ax2 = axes('units', 'centimeters', 'position', ax1.Position);
    maxN = max(dt.NcircM);
    hold on;
    for iA = 1:numel(dt.angleCenters)
        % plot angles and corresponding data
        tmpAngles   = transpose([0, dt.angleCenters(iA) - deg2rad(dt.angularRes / 2):...
            0.001:dt.angleCenters(iA) + deg2rad(dt.angularRes / 2), 0] + dt.rotAngle);
        tmpRadii    = [0; repmat(dt.NcircM(iA), numel(tmpAngles) - 2, 1); 0] ./ maxN;
        [x, y]      = pol2cart(tmpAngles, tmpRadii);
        patch(x, y, [0.7, 0.7, 0.7], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(x, y, '-', ...
            'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    end
    hold off;
    % adjust axes
    set(ax2, ...
        'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
    axis off;
    drawnow;
end

%% Figure 7C

% data
dt  = load(strcat(paths.data, 'Fig_7C_data.mat'));

% figure
figure('units', 'centimeters', 'position', [2, 2, 8, 9]);
axes('units', 'centimeters', 'position', [2.25, 2.8 , 5.6, 6]);
hold on;
for iGroup = 1:numel(dt.groups)
        
    % bar
    bar(iGroup, dt.m{iGroup}, ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.3);
    plot([iGroup, iGroup], [dt.m{iGroup} - dt.ste{iGroup}, dt.m{iGroup} + dt.ste{iGroup}], ...
        'Color', [0, 0, 0], 'LineWidth', 2);
    
    % surrogate distribution
    sp = plot(iGroup + randn(numel(dt.dp{iGroup}), 1) .* 0.05 + 0.2, dt.dp{iGroup}, 'o', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [1, 1, 1], ...
        'MarkerSize', 2);
end
set(gca, ...
    'xlim', [0.4, numel(dt.groups) + 0.6], 'xtick', 1:numel(dt.groups), 'xticklabel', '', ...
    'ylim', [-0.02, 0.03], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
myAxes          = get(gca);
myXTickLabel    = {{'Non-', 'egocentric-', 'bearing'}; {'Egocentric', 'bearing'}};
for iLabel = 1:size(myXTickLabel, 1)
    text(myAxes.XTick(iLabel), min(myAxes.YLim) .* 1.06, myXTickLabel{iLabel}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica', 'FontUnits', 'centimeters', 'FontSize', 0.5);
end
xl = xlabel('Cell type', ...
    'units', 'normalized', 'position', [0.5, -0.35]);
yl = ylabel('Pearson''s \itr');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
lg = legend(sp, 'Surrogates', ...
    'location', 'southwest', 'box', 'off');
lg.ItemTokenSize = [3, 18];
drawnow;

%% Figure 8, B-C, and Figure S8, B-C

% loop through different subpanels
panels  = {'8B'; '8C'; 'S8B'; 'S8C'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], [1, 0, 0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % unsuccessful
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], [0, 0.5, 0], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % successful
    plot(dt.time, dt.m{1}, '-', ...
        'Color', [1, 0, 0], 'LineWidth', 2); % unsuccessful
    plot(dt.time, dt.m{2}, '-', ...
        'Color', [0, 0.5, 0], 'LineWidth', 2); % successful
    set(gca, ...
        'YLim', dt.myYLim);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(dt.timeBorders), max(dt.timeBorders)], 'xtick', min(dt.timeBorders):max(dt.timeBorders), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure 8D and Figure S8D

% loop through different panels
panels  = {'8D'; 'S8D'};
for iPanel = 1:size(panels, 1)

    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], dt.myColor1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % non-spatial cells
    patch([dt.time, fliplr(dt.time)], ...
        [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], dt.myColor2, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % EBCs | spatial cells
    plot(dt.time, dt.m{1}, '-', ...
        'Color', dt.myColor1, 'LineWidth', 2); % non-spatial cells
    plot(dt.time, dt.m{2}, '-', ...
        'Color', dt.myColor2, 'LineWidth', 2); % EBCs | spatial cells
    set(gca, ...
        'YLim', dt.myYLim);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(dt.timeBorders), max(dt.timeBorders)], 'xtick', min(dt.timeBorders):max(dt.timeBorders), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure 8, F-G, and S8, F-G

% loop through different subpanels
panels  = {'8F'; '8G'; 'S8F'; 'S8G'};
for iPanel = 1:size(panels, 1)

    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([dt.time, fliplr(dt.time)], ...
        [dt.mBad - dt.steBad, fliplr(dt.mBad + dt.steBad)], ...
        [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % unsuccessful
    plot(dt.time, dt.mBad, ...
        'Color', [1, 0, 0], 'LineWidth', 2);
    patch([dt.time, fliplr(dt.time)], ...
        [dt.mGood - dt.steGood, fliplr(dt.mGood + dt.steGood)], ...
        [0, 0.5, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % successful
    plot(dt.time, dt.mGood, ...
        'Color', [0, 0.5, 0], 'LineWidth', 2);
    % adjust axis
    set(gca, ...
        'xlim', [min(dt.myXLim), max(dt.myXLim)], 'xtick', min(dt.myXLim):max(dt.myXLim), 'ylim', dt.myYLim, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    % adjust font
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure 8H and Figure S8H

% loop through different panels
panels  = {'8H'; 'S8H'};
for iPanel = 1:size(panels, 1)

    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    patch([dt.thisTime, fliplr(dt.thisTime)], ...
        [dt.m{1} - dt.sem{1}, fliplr(dt.m{1} + dt.sem{1})], dt.myColor1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % non-spatial cells
    patch([dt.thisTime, fliplr(dt.thisTime)], ...
        [dt.m{2} - dt.sem{2}, fliplr(dt.m{2} + dt.sem{2})], dt.myColor2, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3); % EBCs | spatial cells
    plot(dt.thisTime, dt.m{1}, '-', ...
        'Color', dt.myColor1, 'LineWidth', 2); % non-spatial cells
    plot(dt.thisTime, dt.m{2}, '-', ...
        'Color', dt.myColor2, 'LineWidth', 2); % EBCs | spatial cells
    set(gca, ...
        'YLim', dt.myYLim);
    % add x- and y-axis
    myxline(0, '--', [0, 0, 0], 'bottom');
    myyline(0, '--', [0, 0, 0], 'bottom');
    hold off;
    set(gca, ...
        'xlim', [min(dt.myXLim), max(dt.myXLim)], 'xtick', min(dt.myXLim):max(dt.myXLim), ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Time (s)');
    yl = ylabel('Firing rate (Hz)');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S3A and S3B

% loop through different panels
panels  = {'S3A'; 'S3B'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % histogram showing the number of units per wire
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    histogram(dt.numUnitsPerWire, ...
        'FaceColor', [0.7, 0.7, 0.7]);
    set(gca, ...
        'box', 'off', ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    xl = xlabel('Units per wire');
    yl = ylabel('Number of wires');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S3C and S4D

% loop through different panels
panels  = {'S3C'; 'S3D'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % histogram quantifying the ISIs
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    histogram(dt.allPercISILessThan3ms, 0:0.5:20, ...
        'FaceColor', [0.7, 0.7, 0.7]);
    set(gca, ...
        'xlim', [0, 5], 'xtick', 0:1:20, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('% ISI <3 ms');
    yl = ylabel('Number of units');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S3E and S3F

% loop through different panels
panels  = {'S3E'; 'S3F'};
for iPanel = 1:size(panels, 1)
    
    % load data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    histogram(dt.allMeanFR, 0:1:25, ...
        'FaceColor', [0.7, 0.7, 0.7]);
    set(gca, ...
        'xtick', 0:5:25, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Mean FR (Hz)');
    yl = ylabel('Number of units');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S3G and S3H

% loop through different panels
panels  = {'S3G'; 'S3H'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt  = load(strcat(paths.data, 'Fig_', panels{iPanel}, '_data.mat'));
    
    % create figure
    figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    histogram(dt.allPeakSNR, 0:3:33, ...
        'FaceColor', [0.7, 0.7, 0.7]);
    set(gca, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
        'box', 'off');
    xl = xlabel('Peak SNR');
    yl = ylabel('Number of units');
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S4, top

% data
dt  = load(strcat(paths.data, 'Fig_S4Top_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 13.5, 13.5]);

% enhance axes
axes('units', 'normalized', 'position', [0, 0, 1, 1], ...
    'visible', 'off');
text(0.475, 0.015, 'Egocentric bearing', 'units', 'normalized', ...
    'fontunits', 'centimeters', 'fontsize', 0.25);
text(0.01, 0.475, 'Probability', 'units', 'normalized', ...
    'Rotation', 90, ...
    'fontunits', 'centimeters', 'fontsize', 0.25);

% loop trough local reference points
for iLRP = 1:size(dt.locDir.LRPs, 1)
    
    % skip LRPs that were excluded from the analysis
    if dt.locDir.LRPs2Exclude(iLRP, 1) == true
        continue;
    end
            
    % plot average distribution of egocentric direction towards this LRP
    axes('units', 'normalized', 'position', [0.88 * (dt.locDir.LRPs(iLRP, 1) + dt.param.radius) / dt.param.diameter + 0.05, ...
        0.88 * (1 - (dt.locDir.LRPs(iLRP, 2) + dt.param.radius) / dt.param.diameter) + 0.045, 0.065, 0.0625]);
    hold on;
    patch([dt.th{iLRP}, fliplr(dt.th{iLRP})], [dt.m{iLRP} + dt.ste{iLRP}, fliplr(dt.m{iLRP} - dt.ste{iLRP})], [0.7, 0.7, 0.7], ...
        'EdgeColor', 'none');
    plot(dt.th{iLRP}, dt.m{iLRP}, ...
        'Color', [0, 0, 0], 'LineWidth', 1);
    hold off;
    set(gca, ...
        'xlim', [-pi, pi], ...
        'ylim', [0.05, 0.13], ...
        'fontunits', 'centimeters', 'fontsize', 0.17, ...
        'ticklength', [0.05, 0.05]);
    % add xlabel
    if dt.locDir.LRPs2Exclude(iLRP + 1) == true || dt.locDir.LRPs(iLRP, 2) == max(dt.locDir.LRPs(:, 2))
        set(gca, 'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {'B', 'L', 'A', 'R', ''});
    else
        set(gca, 'xtick', []);
    end
    % add ylabel
    if dt.locDir.LRPs(iLRP, 1) == min(dt.locDir.LRPs(:, 1)) || dt.locDir.LRPs2Exclude(iLRP - 12) == true
        set(gca, 'ytick', [0.06, 0.12]);
    else
        set(gca, 'ytick', []);
    end
    % add title
    title([num2str(dt.locDir.LRPs(iLRP, 1)), '/', num2str(dt.locDir.LRPs(iLRP, 2))], ...
        'fontunits', 'centimeters', 'fontsize', 0.17, 'fontweight', 'normal', ...
        'units', 'normalized', 'position', [0.55, 0.85, 0]);    
    drawnow;
end

%% Figure S4, bottom

% data
dt  = load(strcat(paths.data, 'Fig_S4Bottom_data.mat'));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 13.5, 13.5]);

% enhance axes
axes('units', 'normalized', 'position', [0, 0, 1, 1], ...
    'visible', 'off');
text(0.475, 0.01, 'Egocentric bearing', 'units', 'normalized', ...
    'fontunits', 'centimeters', 'fontsize', 0.25);
text(0.01, 0.475, 'Probability', 'units', 'normalized', ...
    'Rotation', 90, ...
    'fontunits', 'centimeters', 'fontsize', 0.25);

% loop trough local reference points
for iLRP = 1:size(dt.locDir.LRPs, 1)
    
    % skip LRPs that were excluded from the analysis
    if dt.locDir.LRPs2Exclude(iLRP, 1) == true
        continue;
    end    
    
    % average distribution of egocentric directions towards this LRP
    axes('units', 'normalized', 'position', [0.88 * (1 - (dt.locDir.LRPs(iLRP, 1) - min(dt.param.XLim)) / dt.param.diameter) + 0.05, ...
        0.88 * (1 - (dt.locDir.LRPs(iLRP, 2) - min(dt.param.ZLim)) / dt.param.diameter) + 0.045, 0.065, 0.0625]);
    hold on;
    patch([dt.ELRPDSamplingPerLRP(iLRP).th, fliplr(dt.ELRPDSamplingPerLRP(iLRP).th)], ...
        [dt.ELRPDSamplingPerLRP(iLRP).m + dt.ELRPDSamplingPerLRP(iLRP).ste, fliplr(dt.ELRPDSamplingPerLRP(iLRP).m - dt.ELRPDSamplingPerLRP(iLRP).ste)], ...
        [0.7, 0.7, 0.7], 'EdgeColor', 'none');
    plot(dt.ELRPDSamplingPerLRP(iLRP).th, dt.ELRPDSamplingPerLRP(iLRP).m, 'Color', [0, 0, 0], 'LineWidth', 1);
    hold off;
    set(gca, ...
        'xlim', [-pi, pi], 'ylim', [0.05, 0.13], ...
        'xdir', 'reverse', ... % because negative angles mean "to the right"
        'fontunits', 'centimeters', 'fontsize', 0.17, ...
        'ticklength', [0.05, 0.05]);
    % add xlabel
    if dt.locDir.LRPs2Exclude(iLRP - 1) == true || dt.locDir.LRPs(iLRP, 2) == max(dt.locDir.LRPs(:, 2))
        set(gca, 'xtick', [-pi, -pi/2, 0, pi/2, pi], 'xticklabel', {'', 'R', 'A', 'L', 'B'}); % because x-axis is reversed
    else
        set(gca, 'xtick', []);
    end
    % add ylabel
    if dt.locDir.LRPs(iLRP, 1) == max(dt.locDir.LRPs(:, 1)) || dt.locDir.LRPs2Exclude(iLRP - dt.locDir.locRes) == true
        set(gca, 'ytick', [0.06, 0.12]);
    else
        set(gca, 'ytick', []);
    end
    % add title
    title([num2str(dt.locDir.LRPs(iLRP, 1)), '/', num2str(dt.locDir.LRPs(iLRP, 2))], ...
        'fontunits', 'centimeters', 'fontsize', 0.17, 'fontweight', 'normal', ...
        'units', 'normalized', 'position', [0.55, 0.85, 0]);
    drawnow;
end

%% Figure S5A (left part)

% data
dt          = load(strcat(paths.data, 'Fig_S5ALeft_data'));
dt          = dt.dt;
dt.visible  = 'on';

% plot figure
LK_PlotFRDir_20200904(dt);
drawnow;

%% Figure S5B

% data
dt  = load(strcat(paths.data, 'Fig_S5B_data.mat'));

% create bar plot for direction cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Direction cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure S5C

% data
dt  = load(strcat(paths.data, 'Fig_S5C_data.mat'));

% create bar plot
figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.75, 2.25, 6, 5.5]);
hold on;
for iBar = 1:numel(dt.m)
    b1 = bar(iBar, dt.m(iBar), ...
        'FaceColor', [0, 0.5, 1]);
    if mod(iBar, 2) == 0
        set(b1, 'FaceColor', [0, 1, 0]);
    elseif mod(iBar, 3) == 0
        set(b1, 'FaceColor', [0.75, 0.75, 0.75]);
    end
    plot([iBar, iBar], [dt.m(iBar) - dt.ste(iBar), dt.m(iBar) + dt.ste(iBar)], ...
        'LineWidth', 2, 'Color', [0 0 0]);
end
% enhance axes
hold off;
set(gca, ...
    'xtick', 1:numel(dt.m), 'xticklabel', '', ...
    'xlim', [0.3, numel(dt.m) + 0.7], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
myAxes          = get(gca);
myXTickLabel    = {'Direction'; {'Egocentric', 'bearing'}};
for iLabel = 1:size(myXTickLabel, 1)
    text(myAxes.XTick(iLabel), -0.02, myXTickLabel{iLabel}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica', 'FontUnits', 'centimeters', 'FontSize', 0.5);
end
xl = xlabel('Cell type', ...
    'Position', [1.5, -0.22]);
yl = ylabel('Vector-field strength (a.u.)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure S6, A-E

% loop through different subpanels
panels  = {'A'; 'B'; 'C'; 'D'; 'E'};
for iPanel = 1:size(panels, 1)
    
    % data
    dt          = load(strcat(paths.data, 'Fig_S6', panels{iPanel}, '_data.mat'));
    dt          = dt.dt;
    dt.visible  = 'on';
    
    % figure for place-like cells
    LK_PlotFRLoc_20200904(dt);
    drawnow;
end

%% Figure S6, F

% data
dt  = load(strcat(paths.data, 'Fig_S6F_data.mat'));

% create bar plot for place-like cells
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
bar(1:numel(dt.uniqueRegions), dt.percPerReg, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, ...
    'box', 'off', 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel('Place-like cells (%)', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;

%% Figure S7, A-D

% loop through the different subpanels
panels  = {'A'; 'B'; 'C'; 'D'};
for iPanel = 1:size(panels, 1)

    % data
    dt  = load(strcat(paths.data, 'Fig_S7', panels{iPanel}, '_data.mat'));
    
    % plot reference points
    figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
    axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
    hold on;
    % boundary and separation between center and periphery
    plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
        'Color', [0, 0, 0], 'LineWidth', 5);
    plot(cos(0:0.001:(2*pi)) .* 2500, sin(0:0.001:(2*pi)) .* 2500, ':', ...
        'Color', [0, 0, 0], 'LineWidth', 2);
    % reference points
    scatter(dt.COMxyCtr(:, 1), dt.COMxyCtr(:, 2), 50, 'o', ...
        'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % center reference points
    scatter(dt.COMxyPeri(:, 1), dt.COMxyPeri(:, 2), 50, 'o', ...
        'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % periphery reference points
    hold off;
    set(gca, ...
        'ydir', 'reverse', ... % negative y-values at top
        'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], 'xtick', [min(dt.locDir.xCenters), max(dt.locDir.xCenters)], ...
        'ylim', [min(dt.locDir.yEdges), max(dt.locDir.yEdges)], 'ytick', [min(dt.locDir.yCenters), max(dt.locDir.yCenters)], ...
        'tickdir', 'out', 'ticklength', [0 0], 'box', 'on');
    xl = xlabel('x (vu)', ...
        'units', 'normalized', 'position', [0.5, -0.025, 0]);
    yl = ylabel('y (vu)', ...
        'units', 'normalized', 'position', [-0.025, 0.5, 0]);
    ytickangle(90);
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    drawnow;
end

%% Figure S7, E

% data
dt  = load(strcat(paths.data, 'Fig_S7E_data.mat'));

% create figure
figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
h1 = histogram(nanmean(dt.Dsurro, 1), ...
    'EdgeColor', 'none', 'FaceColor', [0.7, 0.7, 0.7], ...
    'normalization', 'probability');
plot([nanmean(dt.Demp, 1), nanmean(dt.Demp, 1)], [0, max(h1.Values)], '-', ...
    'Color', [1, 0, 0], 'LineWidth', 2);
xl = xlabel({'Distance between', 'reference points (vu)'});
yl = ylabel('Probability');
set(gca, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
drawnow;
