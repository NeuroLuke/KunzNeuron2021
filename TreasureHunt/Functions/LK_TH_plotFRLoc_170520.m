function f = LK_TH_plotFRLoc_170520(dt)
%
% LK_TH_plotFRLoc_170520 plots firing rate as a function of location.
% 
% Input is a structure with multiple fields.
%
% Output is a figure handle.
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
maxFR = max(dt.FR(:)); % maximum firing rate
imagesc(dt.xCenters, fliplr(dt.zCenters), dt.FR, ...
    'AlphaData', ~isnan(dt.FR));
% colormap
colormap(ax2, 'jet');
% place bins
if isfield(dt, 'locBins')
    locBins                     = dt.locBins .* 0 + maxFR ./ 2;
    locBins(dt.locBins == 1)    = maxFR;
    contour(dt.xCenters, fliplr(dt.zCenters), locBins, 1, ...
        'Color', [0.99, 0.99, 0.99], 'linewidth', 1);
end
% navigation path
if isfield(dt, 'outBeh')
    bMask   = strcmp(dt.outBeh.behInfo(:, 10), 'navigation');
    plot(dt.outBeh.xyz(bMask, 1), dt.outBeh.xyz(bMask, 3), '-', ...
        'Color', [0.7, 0.7, 0.7], ...
        'linewidth', 0.3);
end
hold off;
% flip x- and y-axis
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ...
    'xlim', [min(dt.xCenters), max(dt.xCenters)], 'ylim', [min(dt.zCenters), max(dt.zCenters)]);
axis off;

%----- arena-circle
ax1 = axes('units', 'centimeters', 'position', dt.circleEx);
hold on;
% plot circle
tmpx = -pi:0.001:pi;
patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
    [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [0, 0, 0], ...
    'EdgeColor', 'k');
% xy information
t1 = text(1.1, 0, '(310/360)', ...
    'Rotation', -90);
t2 = text(0, 1.1, '(370/300)');
t3 = text(-1.1, 0, '(430/360)', ...
    'Rotation', 90);
t4 = text(0, -1.1, '(370/420)');
set([t1, t2, t3, t4], ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
    'HorizontalAlignment', 'center');
% indicate maximum firing rate
text(0.175, 0.15, ['max ', num2str(max(max(dt.FR)), '%.2f')], ...
    'units', 'normalized', ...
    'horizontalalignment', 'right', 'verticalalignment', 'top', ...
    'fontunits', 'centimeters', 'fontsize', myFontSize);
hold off;
set(gca, ...
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
    text(0.5, 1.15, ['n=', num2str(sum(dt.numSpikes(dt.bMask4Analysis)))], ...
        'units', 'normalized', ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'HorizontalAlignment', 'center');
end
