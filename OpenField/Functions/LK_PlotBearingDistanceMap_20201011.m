function figHandle = LK_PlotBearingDistanceMap_20201011(cfg)
%
% LK_PlotBearingDistanceMap_20201011 plots a 2D bearing-distance
% firing-rate map.
% 
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2021

%% data preparation

% resize for more accurate plotting of the bearing-distance field
augFac                      = 20;
resBearingBinCenters        = movmean(linspace(min(cfg.bearingBinEdges), max(cfg.bearingBinEdges), augFac * numel(cfg.bearingBinCenters) + 1), 2, 'endpoints', 'discard');
resDistanceBinCenters       = movmean(linspace(min(cfg.distanceBinEdges), max(cfg.distanceBinEdges), augFac * numel(cfg.distanceBinCenters) + 1), 2, 'endpoints', 'discard');
resBearingDistanceField     = myresizem(cfg.bearingDistanceField, augFac);

% indices of distance edges surrounding distance bins with firing rates
meanDistanceTuning  = nanmean(cfg.firingRate, 2);
idxValidDistance    = find(~isnan(meanDistanceTuning), 1, 'first'):find(~isnan(meanDistanceTuning), 1, 'last') + 1;

%% 2D bearing-distance plot

% create figure
figHandle = figure('units', 'centimeters', 'position', [2, 2, 9, 7], ...
    'visible', cfg.visible);

% axes for 2D bearing-distance map
ax1 = axes('units', 'centimeters', 'position', [1.5, 1.6, 4, 4]);
hold on;
% plot firing-rate map and the bearing-distance field
imagesc(cfg.bearingBinCenters, cfg.distanceBinCenters, cfg.firingRate, ...
    'alphadata', ~isnan(cfg.firingRate));
contour(resBearingBinCenters, resDistanceBinCenters, resBearingDistanceField, 1, ...
    'Color', [0, 0, 0], 'LineWidth', 1);
hold off;
colormap jet;
cb = colorbar;
caxis([min(cfg.firingRate(:)), quantile(cfg.firingRate(:), 0.995)]);
ylabel(cb, 'Hz', 'Rotation', 0, ...
    'units', 'normalized', 'position', [1.75, 0.5], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
cb.Units = 'centimeters';
cb.Position = [6, 1, 0.25, 1.5];
cb.Limits = round(cb.Limits, 1);
cb.Ticks = round(cb.Limits, 1);
cb.Label.FontUnits = 'centimeters';
cb.Label.FontSize = 0.4;
xl = xlabel('Bearing');
yl = ylabel('Distance (vu)', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
myT = title(cfg.figTitle, ...
    'units', 'normalized', 'position', [0.5, 1.02, 0]);
set(gca, ...
    'xlim', [min(cfg.bearingBinEdges), max(cfg.bearingBinEdges)], 'xtick', linspace(min(cfg.bearingBinEdges), max(cfg.bearingBinEdges), 5), ...
    'xticklabel', {'B', 'L', 'A', 'R', ''}, ...
    'ylim', [0, max(cfg.distanceBinEdges(idxValidDistance))], 'ytick', [min(cfg.distanceBinEdges), max(cfg.distanceBinEdges(idxValidDistance))], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl, myT], ...
    'Fontunits', 'centimeters', 'Fontsize', 0.5, 'FontWeight', 'normal');

%% waveform

% axes for waveform
ax2 = axes('units', 'centimeters', 'position', [6.7, 4, 2, 2]);
hold on;
% spike density plot (Reber et al., 2019)
spikeTime   = (0:size(cfg.thisSpike, 2) - 1) ./ cfg.sr .* 1000; % (ms)
outPlot     = LK_DensityPlot(spikeTime, cfg.thisSpike);
hold off;
% enhance axes
set(gca, ...
    'xlim', round([min(spikeTime), max(spikeTime)]), 'xtick', round([min(spikeTime), max(spikeTime)]), ...
    'ylim', [outPlot.lbound, outPlot.ubound], 'ytick', [outPlot.lbound, outPlot.ubound], ...
    'ticklength', [0, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.4);
colormap(ax2, 'parula');
% enhance axes
xl = xlabel('ms', ...
    'units', 'normalized', 'position', [0.5, -0.1, 0]);
yl = ylabel('\muV', ...
    'units', 'normalized', 'position', [-0.12, 0.5, 0]);
tx = text(0.5, 1.15, ['n=', num2str(cfg.nspk)], ...
    'units', 'normalized', ...
    'HorizontalAlignment', 'center');
set([xl, yl, tx], ...
    'FontUnits', 'centimeters', 'FontSize', 0.4)
