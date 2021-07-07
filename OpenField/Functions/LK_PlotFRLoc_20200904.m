function f = LK_PlotFRLoc_20200904(dt)
%
% LK_PlotFRLoc_20200904 plots firing rate as a function of location.
%
% Input is a structure with multiple fields including
%   FR              --> firing rate
%   xCenters        --> x centers
%   yCenters        --> y centers
%   locBins         --> location bins
%   ...
%
% Output is a handle to the figure.
%
% Lukas Kunz, 2021

% font size
myFontSize  = 0.4;

%% create figure

% basics
f = figure('units', 'centimeters', 'position', dt.figEx);
if strcmp(dt.visible, 'off')
    set(f, 'visible', 'off');
end
axis off;

%----- firing rate
ax2 = axes('units', 'centimeters', 'position', dt.circleEx);
hold on;
% firing rate map
maxFR = max(max(dt.FR)); % maximum firing rate
imagesc(dt.xCenters, fliplr(dt.yCenters), dt.FR, ...
    'AlphaData', ~isnan(dt.FR)); % low y-values are at the bottom
% colormap
colormap(ax2, 'jet');
% place bins
if isfield(dt, 'locBins')
    locBins                     = dt.locBins .* 0 + maxFR ./ 2;
    locBins(dt.locBins == 1)    = maxFR;
    contour(dt.xCenters, fliplr(dt.yCenters), locBins, 1, ...
        'Color', [0.99, 0.99, 0.99], 'linewidth', 1);
end
% plot navigation path
if isfield(dt, 'path')
    plot(dt.path(:, 1), dt.path(:, 2), '-', ...
        'Color', [0.7, 0.7, 0.7], 'linewidth', 0.3);
end
hold off;
% flip y-axis
set(gca, ...
    'ydir', 'reverse', ...
    'xlim', [-1.3, 1.3] .* max(dt.xCenters), 'ylim', [-1.3, 1.3] .* max(dt.yCenters));
axis off;

%----- boundary
ax1 = axes('units', 'centimeters', 'position', dt.circleEx);
hold on;
% plot circle
tmpx = -pi:0.001:pi;
patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
    [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [0, 0, 0], ...
    'EdgeColor', 'k');
% plot xy-value text
t1 = text(1.1, 0, '(5000/0)', ...
    'Rotation', -90);
t2 = text(0, 1.1, '(0/5000)');
t3 = text(-1.1, 0, '(-5000/0)', ...
    'Rotation', 90);
t4 = text(0, -1.1, '(0/-5000)');
set([t1, t2, t3, t4], ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
    'HorizontalAlignment', 'center');
% indicate maximum firing rate
text(0.175, 0.15, [num2str(max(max(dt.FR)), '%.1f'), ' Hz'], ...
    'units', 'normalized', ...
    'horizontalalignment', 'right', 'verticalalignment', 'top', ...
    'fontunits', 'centimeters', 'fontsize', myFontSize);
hold off;
set(gca, ...
    'ydir', 'reverse', ... % negative y-values are pointing north
    'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
axis off;
title(dt.figTitle, ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.1, 0]);

%----- waveform
if isfield(dt, 'waveEx')
    
    ax3 = axes('units', 'centimeters', 'position', dt.waveEx);    
    hold on;
    
    % spike density plot (Reber et al, 2019)
    spikeTime   = (0:size(dt.thisSpike, 2) - 1) ./ dt.t.par.sr .* 1000;
    outPlot     = LK_DensityPlot(spikeTime, dt.thisSpike);
    hold off;
    
    % enhance axes
    set(gca, ...
        'xlim', round([min(spikeTime), max(spikeTime)]), 'xtick', round([min(spikeTime), max(spikeTime)]), ...
        'ylim', [outPlot.lbound, outPlot.ubound], 'ytick', [outPlot.lbound, outPlot.ubound], ...
        'ticklength', [0, 0], ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize);
    xlabel('ms', ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'units', 'normalized', 'position', [0.5, -0.1, 0]);
    ylabel('\muV', ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'units', 'normalized', 'position', [-0.12, 0.5, 0]);
    text(0.5, 1.15, ['n=', num2str(dt.nspk)], ...
        'units', 'normalized', ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'HorizontalAlignment', 'center');
end
