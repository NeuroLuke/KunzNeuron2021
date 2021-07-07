function f = LK_PlotFRObj_20210515(dt)
%
% LK_PlotFRObj_20210515 plots firing rate as a function of object.
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

%----- firing-rate
axes('units', 'centimeters', 'position', dt.boxEx);
hold on;
for iB = 1:numel(dt.FR)
    % mean
    myb = bar(iB, dt.FR(iB), ...
        'FaceColor', [1, 1, 1]);
    if dt.objBins(iB) == true
        set(myb, 'FaceColor', rgb('orange')); % preferred object
    end
    % SEM
    plot([iB, iB], [dt.FR(iB) - dt.SEM(iB), dt.FR(iB) + dt.SEM(iB)], ...
        'Color', [0, 0, 0], 'LineWidth', 2);
end
hold off;
maxY = ceil(max(dt.FR) .* 1.1);
set(gca, ...
    'xlim', [0, numel(dt.FR) + 1], 'xtick', 1:numel(dt.FR), ...
    'ylim', [0, maxY], 'ytick', [0, maxY], ...
    'fontunits', 'centimeters', 'fontsize', myFontSize, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xlabel('Object', ...
    'fontunits', 'centimeters', 'fontsize', myFontSize);
ylabel('FR (Hz)', ...
    'fontunits', 'centimeters', 'fontsize', myFontSize, ...
    'units', 'normalized', 'position', [-0.04, 0.5, 0]);
title(dt.figTitle, ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.1, 0]);

%----- waveform
if isfield(dt, 'waveEx')
    axes('units', 'centimeters', 'position', dt.waveEx);
    hold on;
    
    % density plot (Reber et al., 2019)
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

%--- object-locations
axes('units', 'centimeters', 'position', [0.5, 0.5, 2.5, 2.5]);
hold on;
% plot object locations
for iObj = min(dt.object.idx, 1):max(dt.object.idx)
    
    % this object location
    thisObjXY   = unique(dt.trials(dt.trials(:, 2) == iObj, 9:10), 'rows', 'stable');
    thisObjXY   = thisObjXY(~isnan(thisObjXY(:, 1)), :);
    myp         = plot(thisObjXY(1, 1), thisObjXY(1, 2), 'ko');
    
    % preferred objects
    if dt.objBins(iObj + 1) == 1
        set(myp, 'MarkerFaceColor', rgb('orange')); % preferred object
    end
end
% plot boundary
plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, 'k-', ...
    'LineWidth', 2);
t1 = text(0, 6000, 'x');
t2 = text(-6000, 0, 'y');
t3 = text(0, -7750, {'Object', 'locations'});
set([t1, t2, t3], ...
    'fontunits', 'centimeters', 'fontsize', myFontSize, ...
    'HorizontalAlignment', 'center');
hold off;
set(gca, ...
    'ydir', 'reverse', ...
    'xlim', [-5000, 5000], 'ylim', [-5000, 5000]);
axis off square equal;
