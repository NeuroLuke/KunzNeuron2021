function f = LK_PlotFRDir_20200904(dt)
%
% LK_PlotFRDir_20200904 plots firing rate as a function of direction.
%
% Input is a structure with multiple fields including
%   FR              --> firing rate
%   angleCenters    --> angle centers
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

%----- angle-circle
ax1 = axes('units', 'centimeters', 'position', dt.circleEx);
hold on;
% plot angle lines
for iAng = 0:30:165
    plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
        'Color', rgb('gray'), 'LineWidth', 1);
end
% plot circle
tmpx = -pi:0.001:pi;
tmpc = linspace(0, 1, length(tmpx));
patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
    [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
    'EdgeColor', 'interp');
% plot angle text
t1 = text(1.05, 0, dt.angLabels{1}, ...
    'HorizontalAlignment', 'left');
t2 = text(0, 1.125, dt.angLabels{2}, ...
    'HorizontalAlignment', 'center');
t3 = text(-1.05, 0, dt.angLabels{3}, ...
    'HorizontalAlignment', 'right');
t4 = text(0, -1.125, dt.angLabels{4}, ...
    'HorizontalAlignment', 'center');
% indicate maximum firing rate
t5 = text(0.175, 0.15, [num2str(max(dt.FR), '%.1f'), ' Hz'], ...
    'units', 'normalized', ...
    'horizontalalignment', 'right', 'verticalalignment', 'top');
set([t1, t2, t3, t4, t5], ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize);
hold off;
set(gca, ...
    'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
% flip y-axis because the arena is programmed such that negative yaw-values
% point north
if dt.bReverseY == true
    set(gca, ...
        'ydir', 'reverse');
end
myCM = circshift(hsv, size(hsv, 1) / 2); % so that 0° is red
colormap(ax1, myCM);
axis off;

%----- firing rate
ax2 = axes('units', 'centimeters', 'position', ax1.Position);
maxFR = max(dt.FR);
hold on;
for iA = 1:size(dt.angleCenters, 1)
    thisAngles  = transpose([0, dt.angleCenters(iA) - deg2rad(dt.angularRes / 2):0.001:dt.angleCenters(iA) + deg2rad(dt.angularRes / 2), 0]);
    thisRadii   = [0; repmat(dt.FR(iA), numel(thisAngles) - 2, 1); 0] ./ maxFR; % normalize in relation to maximum FR
    [x, y]      = pol2cart(thisAngles, thisRadii);
    patch(x, y, [0.7, 0.7, 0.7], ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    p1 = plot(x, y, '-', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    if dt.angularRes < 15
        set(p1, 'LineWidth', 0.5);
    end
end
% circular mean
circM               = circ_mean(dt.angleCenters, dt.FR);
[circMx, circMy]    = pol2cart(circM, 1);
plot([0, circMx], [0, circMy], '-', ...
    'linewidth', 3, 'color', [0, 0, 0]);
hold off;
set(gca, ...
    'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
if dt.bReverseY == true
    set(gca, ...
        'ydir', 'reverse'); % flip y-axis
end
axis off;
title(dt.figTitle, ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.1, 0]);

%----- waveform
if isfield(dt, 'waveEx')
    ax3 = axes('units', 'centimeters', 'position', dt.waveEx);
    hold on;
    
    % spike density plot (Reber et al., 2019)
    spikeTime   = (0:size(dt.thisSpike, 2) - 1) ./ dt.t.par.sr .* 1000; % (ms)
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
    % indicate number of spikes used for the analysis
    text(0.5, 1.15, ['n=', num2str(dt.nspk)], ...
        'units', 'normalized', ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'HorizontalAlignment', 'center');
end
