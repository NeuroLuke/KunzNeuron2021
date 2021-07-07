function f = LK_TH_PlotELRPDModulation_20210515(dt)
%
% LK_TH_PlotELRPDModulation_20210515 plots a depiction of FR modulation by
% egocentric local reference point direction.
%
% Input is a structure with fields
% 	visible            	--> visibility of figure ('on' or 'off')
%   param               --> general settings
%  	locDir              --> settings for ELRPD
% 	locdir_CBPT       	--> output from cluster-based permutation test
%  	all_locdir_corrFR   --> firing rate for ELRPD
% 	direc               --> settings for direction
% 	dir_corrFR_perLRP   --> location-specific directional firing rate
% 	thisSpike           --> spike waveforms
% 	t.par.sr            --> sampling rate
% 	nspk                --> number of spikes
% 	figTitle            --> figure title
%
% Output is a figure handle.
%
% Lukas Kunz, 2021

% font size
myFontSize  = 0.4;

% rotate the angles (so that "ahead" = 0° is at the top)
rotAngle    = deg2rad(90); % positive values rotate counterclockwise

%% create figure

% basis
f = figure('units', 'centimeters', 'position', [2, 5, 22.8, 7.8]);
if strcmp(dt.visible, 'off')
    set(f, 'visible', 'off');
end

% arena boundary
ax1 = axes('units', 'centimeters', 'position', [0.8, 1, 5.5, 5.5]);
plot(dt.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* dt.param.arenaRadius, dt.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* dt.param.arenaRadius, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% flip x- and y-axis
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], 'xtick', [min(dt.locDir.xCenters), max(dt.locDir.xCenters)], ...
    'ylim', [min(dt.locDir.zEdges), max(dt.locDir.zEdges)], 'ytick', [min(dt.locDir.zCenters), max(dt.locDir.zCenters)], ...
    'ticklength', [0 0], ...
    'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('z (vu)', ...
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
myColors   	= flipud(circshift(hsv, size(hsv, 1) / 2)); % so that red means "ahead" and "right" is associated with green
potAngles  	= linspace(-pi, pi, size(myColors, 1) + 1);
for iLRP = 1:size(dt.locDir.LRPs, 1)
    
    % if this LRP is not included in the analyses
    if dt.locDir.LRPs2Exclude(iLRP) == true
        continue;
    end
    
    % circular mean of tuning curve at this local reference point
    bIsNan  = isnan(dt.all_locdir_corrFR(iLRP, :));
    if sum(~bIsNan) >= 2
        
        % preferred ELRPD for this LRP
        prefAng     = circ_mean(dt.locDir.angleCenters(~bIsNan)', dt.all_locdir_corrFR(iLRP, ~bIsNan)'); % circular mean
        prefAngBin  = discretize(prefAng, potAngles);
        p1 = plot(dt.locDir.LRPs(iLRP, 1), dt.locDir.LRPs(iLRP, 2), 'o', ...
            'MarkerSize', 4, ...
            'MarkerFaceColor', [0.7, 0.7, 0.7], ...
            'MarkerEdgeColor', 'none');
        
        % if this LRP is significant
        if dt.locdir_CBPT.largestCluster(iLRP, 1) > 0.95
            set(p1, ...
                'MarkerSize', 10, 'MarkerFaceColor', myColors(prefAngBin, :));
        end
    end
end
% indicate LRP for which the ELRPD tuning is shown in detail
if ~isempty(dt.locdir_CBPT.COMxzClosestLRP)
    plot(dt.locdir_CBPT.COMxzClosestLRP(1, 1), dt.locdir_CBPT.COMxzClosestLRP(1, 2), 'o', ...
        'Color', [0, 0, 0], ...
        'MarkerSize', 10, 'LineWidth', 2);
end
hold off;
% flip x- and y-axis
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], ...
    'ylim', [min(dt.locDir.zEdges), max(dt.locDir.zEdges)]);

%--- relationship between FR and ELRPD for COM
if ~isempty(dt.locdir_CBPT.COMxzClosestLRP)
    
    % index of LRP closest to the reference point
    tmpMaxIdx   = find(all(dt.locdir_CBPT.COMxzClosestLRP == dt.locDir.LRPs, 2));
    
    %----- angle-circle
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
    % angle information
    t1 = text(1.1125, 0, 'R'); % "east" of figure
    t2 = text(0, 1.125, 'A'); % "north" of figure
    t3 = text(-1.1, 0, 'L');
    t4 = text(0, -1.125, 'B');
    set([t1, t2, t3, t4], ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
        'HorizontalAlignment', 'center');
    % indicate maximum firing rate
    text(0.825, 0.15, [num2str(max(dt.all_locdir_corrFR(tmpMaxIdx, :)), '%.1f'), ' Hz'], ...
        'units', 'normalized', ...
        'horizontalalignment', 'left', 'verticalalignment', 'top', ...
        'fontunits', 'centimeters', 'fontsize', myFontSize)
    hold off;
    set(gca, ...
        'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
    colormap(ax4, myColors);
    axis off;
    
    %----- firing rate
    ax5 = axes('units', 'centimeters', 'position', ax4.Position);
    thisFR  = transpose(dt.all_locdir_corrFR(tmpMaxIdx, :));
    maxFR  	= max(thisFR);
    hold on;
    for iA = 1:numel(dt.locDir.angleCenters)
        thisAngles  = transpose([0, dt.locDir.angleCenters(iA) - deg2rad(dt.locDir.angularRes / 2):...
            0.001:dt.locDir.angleCenters(iA) + deg2rad(dt.locDir.angularRes / 2), 0]);
        thisRadii  	= [0; repmat(thisFR(iA), numel(thisAngles) - 2, 1); 0] ./ maxFR; % normalize in relation to maximum FR
        [x, y]   	= pol2cart(thisAngles + rotAngle, thisRadii); % rotate angles by +90°
        patch(x, y, [0.7, 0.7, 0.7], ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        plot(x, y, '-', ...
            'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    end
    % circular mean
    circM               = circ_mean(dt.locDir.angleCenters', thisFR);
    [circMx, circMy]    = pol2cart(circM + rotAngle, 1); % rotate angles by +90°
    plot([0, circMx], [0, circMy], '-', ...
        'linewidth', 3, 'color', [0, 0, 0]);    
    hold off;
    set(gca, ...
        'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
    axis off;
    title({'Tuning curve', 'for reference point'}, ...
        'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
        'Units', 'normalized', 'Position', [0.5, 1.075, 0]);
end

%--- waveform
ax6 = axes('units', 'centimeters', 'position', [13.4, 5, 2, 2]);
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
text(0.5, 1.15, ['n=', num2str(dt.nspk)], ...
    'units', 'normalized', ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
    'HorizontalAlignment', 'center');

%--- quiver plot
% for each LRP, extract the preferred allocentric heading direction ("PD")
meanPDPerLRP    = nan(size(dt.dir_corrFR_perLRP, 1), 1);
for iLRP = 1:size(dt.dir_corrFR_perLRP, 1)
    bNaN    = isnan(dt.dir_corrFR_perLRP(iLRP, :));
    if sum(~bNaN) >= 2
        meanPDPerLRP(iLRP)  = circ_mean(dt.direc.angleCenters(~bNaN)', transpose(dt.dir_corrFR_perLRP(iLRP, ~bNaN)));
    end
end
% plot boundary
ax7 = axes('units', 'centimeters', 'position', [16.9, 1, 5.5, 5.5]);
plot(dt.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* dt.param.arenaRadius, dt.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* dt.param.arenaRadius, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
myT = title({'Preferred heading', 'directions'}, ...
    'units', 'normalized', 'position', [0.5, 1.02, 0]);
% flip x- and y-axis
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], 'xtick', [min(dt.locDir.xCenters), max(dt.locDir.xCenters)], ...
    'ylim', [min(dt.locDir.zEdges), max(dt.locDir.zEdges)], 'ytick', [min(dt.locDir.zCenters), max(dt.locDir.zCenters)], ...
    'ticklength', [0, 0], ...
    'box', 'on');
ytickangle(ax7, 90); % rotate y-tick-labels
set([gca, xl, yl, myT], ...
    'fontunits', 'centimeters', 'fontsize', myFontSize, 'fontweight', 'normal');

% plot preferred heading-direction modulation for each LRP
ax8 = axes('units', 'centimeters', 'position', ax7.Position);
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ...
    'xlim', [min(dt.locDir.xEdges), max(dt.locDir.xEdges)], ...
    'ylim', [min(dt.locDir.zEdges), max(dt.locDir.zEdges)]);
hold on;
% convert preferred directions into Cartesian coordinates
[qu, qv]    = pol2cart(meanPDPerLRP, ones(size(meanPDPerLRP))); % set radius to one for all LRPs
lineLength  = 7.5;
for iLRP = 1:size(dt.locDir.LRPs, 1)
    
    % if you do not want to use this LRP or if there is no data
    if dt.locDir.LRPs2Exclude(iLRP) == true || isnan(qu(iLRP, 1))
        continue;
    end
    
    % plot directional tuning for this LRP
    arrow([dt.locDir.LRPs(iLRP, 1), dt.locDir.LRPs(iLRP, 2)], [dt.locDir.LRPs(iLRP, 1) + qu(iLRP, 1) .* lineLength, dt.locDir.LRPs(iLRP, 2) + qv(iLRP, 1) .* lineLength], ...
        7, 'BaseAngle', 40, 'TipAngle', 25, 'width', 1.1, ...
        'Color', [0.7, 0.7, 0.7]);
end
hold off;
axis off;
