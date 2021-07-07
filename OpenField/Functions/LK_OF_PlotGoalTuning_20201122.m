function f = LK_OF_PlotGoalTuning_20201122(dt)
%
% LK_OF_PlotGoalTuning_20201122 plots a depiction of FR modulation by
% egocentric local reference point direction.
%
% Input is a structure with multiple fields.
% 
% Output is a figure handle.
%
% Lukas Kunz, 2021

% font size
myFontSize  = 0.4;

% determine whether you want to rotate the angles (so that "ahead" = 0° is
% at the top)
rotAngle    = deg2rad(90); % positive values rotate counter-clockwise

%% create figure

% identify goal with strongest egocentric tuning
[~, idxMaxTuning]   = max(dt.locdir_Frank);

% basis
f = figure('units', 'centimeters', 'position', [2, 5, 15.5, 8]);
if strcmp(dt.visible, 'off')
    set(f, 'visible', 'off');
end

% arena boundary
ax1 = axes('units', 'centimeters', 'position', [0.8, 1, 5.5, 5.5]);
plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% flip y-axis because, from bird's eyeview, negative y-values are at top
set(gca, ...
    'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], 'xtick', dt.locDir.xCenters([1, end]), ...
    'ylim', [min(dt.locDir.yEdges), max(dt.locDir.yEdges)], 'ytick', dt.locDir.yCenters([1, end]), ...
    'ticklength', [0, 0], 'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
myT = title(dt.figTitle, ...
    'units', 'normalized', 'position', [0.5, 1.02, 0]);
set([gca, xl, yl, myT], ...
    'fontunits', 'centimeters', 'fontsize', myFontSize', 'fontweight', 'normal');
ytickangle(ax1, 90); % rotate y-tick-labels

% F-rank map together with information about preferred ELRPD
ax3 = axes('units', 'centimeters', 'position', ax1.Position, ...
    'visible', 'off');
hold on;
myColors  	= circshift(hsv, size(hsv, 1) / 2); % so that red means ahead and +90° is associated with green
potAngles 	= linspace(-pi, pi, size(myColors, 1) + 1);
for iLRP = 1:size(dt.locDir.LRPs, 1)
    
    % circular mean of tuning curve at this local reference point
    bIsNan  = isnan(dt.locdir_corrFR(iLRP, :));
    if sum(~bIsNan) >= 2
        
        % preferred ELRPD for this LRP
        prefAng     = circ_mean(dt.locDir.angleCenters(~bIsNan), dt.locdir_corrFR(iLRP, ~bIsNan)'); % circular mean
        prefAngBin  = discretize(prefAng, potAngles);
        p1 = plot(dt.locDir.LRPs(iLRP, 1), dt.locDir.LRPs(iLRP, 2), 'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', myColors(prefAngBin, :), ...
            'MarkerEdgeColor', 'none');
        
        % if this LRP is a response location
        if iLRP >= 9
            set(p1, 'Marker', 'p', 'MarkerSize', 8);
        end
        
        % if this LRP is significant
        if dt.locdir_Frank(iLRP) > 0.95
            set(p1, 'MarkerSize', 10); % double the size
            if iLRP >= 9
                set(p1, 'MarkerSize', 16); % double the size
            end
        end
    end
end
% indicate LRP for which the ELRPD tuning is depicted in detail
plot(dt.locDir.LRPs(idxMaxTuning, 1), dt.locDir.LRPs(idxMaxTuning, 2), 'o', ...
    'Color', [0, 0, 0], ...
    'MarkerSize', 10, 'LineWidth', 2);
hold off;
% flip y-axis because, from bird's eyeview, negative y-values are at top
set(gca, ...
    'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], ...
    'ylim', [min(dt.locDir.yEdges), max(dt.locDir.yEdges)]);

%--- relationship between FR and ELRPD for most strongly modulated goal

%- angle-circle
ax4 = axes('units', 'centimeters', 'position', [7, 0.75, 5.5, 5.5]);
hold on;
% plot angle lines
for iAng = 0:30:165
    plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
        'Color', rgb('gray'), 'LineWidth', 1);
end
% plot circle
tmpx = (-pi:0.001:pi) + rotAngle; % rotate angles so that ahead is at the top
tmpc = linspace(0, 1, length(tmpx));
patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
    [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
    'EdgeColor', 'interp');
% plot angle text
t1 = text(1.15, 0, 'L'); % L because x-axis will be reversed
t2 = text(0, 1.125, 'A');
t3 = text(-1.15, 0, 'R'); % R because x-axis will be reversed
t4 = text(0, -1.125, 'B');
set([t1, t2, t3, t4], ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
    'HorizontalAlignment', 'center');
% indicate maximum firing rate
text(0.825, 0.15, [num2str(max(dt.locdir_corrFR(idxMaxTuning, :)), '%.1f'), ' Hz'], ...
    'units', 'normalized', ...
    'horizontalalignment', 'left', 'verticalalignment', 'top', ...
    'fontunits', 'centimeters', 'fontsize', myFontSize)
hold off;
% cave: reverse x-axis so that +90° ("right") is right
set(gca, ...
    'xdir', 'reverse', ...
    'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
colormap(ax4, myColors);
axis off;

%- firing rate
ax5 = axes('units', 'centimeters', 'position', ax4.Position);
thisFR  = transpose(dt.locdir_corrFR(idxMaxTuning, :));
maxFR  	= max(thisFR);
hold on;
for iA = 1:size(dt.locDir.angleCenters, 1)
    thisAngles  = transpose([0, dt.locDir.angleCenters(iA) - deg2rad(dt.locDir.angularRes / 2):...
        0.001:dt.locDir.angleCenters(iA) + deg2rad(dt.locDir.angularRes / 2), 0]);
    thisRadii   = [0; repmat(thisFR(iA), numel(thisAngles) - 2, 1); 0] ./ maxFR; % normalize in relation to maximum FR
    [x, y]     	= pol2cart(thisAngles + rotAngle, thisRadii); % rotate angles by +90°
    patch(x, y, [0.7, 0.7, 0.7], ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(x, y, '-', ...
        'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
end
% circular mean
circM               = circ_mean(dt.locDir.angleCenters, thisFR);
[circMx, circMy]    = pol2cart(circM + rotAngle, 1); % rotate angles by +90°
plot([0, circMx], [0, circMy], '-', ...
    'linewidth', 3, 'color', [0, 0, 0]);
hold off;
% cave: reverse x-axis so that +90° is right
set(gca, ...
    'xdir', 'reverse', ...
    'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
axis off;
title({'Tuning curve', 'for strongest goal'}, ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.075, 0]);

%--- waveform
ax6 = axes('units', 'centimeters', 'position', [13.25, 5.3, 2, 2]);
hold on;
% spike density plot (Reber et al., 2019)
spikeTime   = ((0:size(dt.thisSpike, 2) - 1) ./ dt.t.par.sr) .* 1000; % (msec)
outPlot     = LK_DensityPlot(spikeTime, dt.thisSpike);
hold off;
set(gca, ...
    'xlim', round([min(spikeTime), max(spikeTime)]), 'xtick', round([min(spikeTime), max(spikeTime)]), ...
    'ylim', [outPlot.lbound, outPlot.ubound], 'ytick', [outPlot.lbound, outPlot.ubound], ...
    'ticklength', [0, 0], ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize);
colormap(ax6, 'parula');
% enhance axes
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

