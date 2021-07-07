function [f, linFit] = LK_PlotMemoryTuning_20210515(dt)
%
% LK_PlotMemoryTuning_20210515 plots the relationship between memory
% performance and firing rates.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle and the linear fit.
%
% Lukas Kunz, 2021

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8.5, 6.5], ...
    'visible', dt.visible);

%----- memory effect on firing rate
ax1 = axes('units', 'centimeters', 'position', [4.7, 1.4, 3.5, 3.8]);
hold on;
plot(dt.memorLevels, dt.empValsCorr, 'o', ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerSize', 5); % data
linFit = polyfit(dt.memorLevels, dt.empValsCorr, 1); % linear fit
tmpAx = get(gca);
plot(tmpAx.XLim, tmpAx.XLim .* linFit(1) + linFit(2), '-', ...
    'Color', [1, 0, 0], 'LineWidth', 2);
hold off;
set(gca, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
if strcmp(dt.memor.bRankData, 'yes')
    xl = xlabel('Performance (ranked)');
else
    xl = xlabel('Performance');
end
yl = ylabel('FR residuals (Hz)');
tl = title(dt.figTitle, ...
    'units', 'normalized', 'position', [0.5, 1.025]);
set([gca, xl, yl, tl], ...
    'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');

%----- waveform
ax2 = axes('units', 'centimeters', 'position', [1, 4, 2, 2]);
hold on;
% spike density plot (Reber et al., 2019)
spikeTime   = (0:size(dt.thisSpike, 2) - 1) ./ dt.t.par.sr .* 1000;
outPlot     = LK_DensityPlot(spikeTime, dt.thisSpike);
hold off;
% enhance axes
set(gca, ...
    'xlim', round([min(spikeTime), max(spikeTime)]), 'xtick', round([min(spikeTime), max(spikeTime)]), ...
    'ylim', [outPlot.lbound, outPlot.ubound], 'ytick', [outPlot.lbound, outPlot.ubound], ...
    'ticklength', [0, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.4);
xlabel('ms', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'units', 'normalized', 'position', [0.5, -0.1, 0]);
ylabel('\muV', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'units', 'normalized', 'position', [-0.12, 0.5, 0]);
text(0.5, 1.15, ['n=', num2str(dt.nspk)], ...
    'units', 'normalized', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'HorizontalAlignment', 'center');
