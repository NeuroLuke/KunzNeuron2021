function [f, P]	= LK_PlotLinearDistModulation_20201011(cfg)
%
% LK_PlotLinearDistModulation_20201011 plots the linear distance modulation
% of firing rates.
%
% Input is a structure with fields
%   visible         --> 'on' or 'off'
%   FR              --> firing rates
%   distnc          --> distance information
%   locdir_COMxy    --> reference point
%   prefELRPD       --> preferred egocentric bearing
%   thisSpike       --> spike waveforms
%   sr              --> sampling rate
%   nspk            --> number of spikes
%   figTitle        --> figure title
%
% Output is a figure handle and the parameters for the linear fit.
%
% Lukas Kunz, 2021

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 6.5], ...
    'visible', cfg.visible);

%----- distance effect on firing rate
ax1 = axes('units', 'centimeters', 'position', [4.3, 1.5, 3.25, 3.7]);
hold on;
bValid = ~isnan(cfg.FR);
bar(cfg.distnc.binCenters(bValid), cfg.FR(bValid), ...
    'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
P = polyfit(cfg.distnc.binCenters(bValid), cfg.FR(bValid), 1);
thisX = [-250, max(cfg.distnc.binCenters(bValid)) + 500];
plot(thisX, thisX .* P(1) + P(2), '-', ...
    'Color', [1, 0, 0], 'LineWidth', 2);
hold off;
xl = xlabel('Distance (vu)');
yl = ylabel('FR (Hz)', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
tl = title(cfg.figTitle, ...
    'units', 'normalized', 'position', [0.5, 1.025]);
set([gca, xl, yl, tl], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5, 'FontWeight', 'normal');
tmpAx = get(gca);
set(gca, ...
    'xlim', thisX, 'xtick', [min(cfg.distnc.binCenters), max(cfg.distnc.binCenters(bValid))], ...
    'ylim', [min(tmpAx.YLim), max(tmpAx.YLim)], 'ytick', [myceil(min(tmpAx.YLim), 1), myfloor(max(tmpAx.YLim), 1)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');

%----- reference point and associated preferred ELRPD
ax2 = axes('units', 'centimeters', 'position', [0.9, 0.5, 2, 2]);
hold on;
% plot circles for distance-bin centers
for iC = 1:numel(cfg.distnc.binCenters)
    if ~bValid(iC)
        continue;
    end
    tmpX                = cfg.locdir_COMxy(1) + cos(0:0.001:(2*pi)) .* cfg.distnc.binCenters(iC);
    tmpY                = cfg.locdir_COMxy(2) + sin(0:0.001:(2*pi)) .* cfg.distnc.binCenters(iC);
    [tmpTheta, tmpRho]  = cart2pol(tmpX, tmpY);
    tmpX(tmpRho > 5000) = cos(tmpTheta(tmpRho > 5000)) .* 5000;
    tmpY(tmpRho > 5000) = sin(tmpTheta(tmpRho > 5000)) .* 5000;
    plot(tmpX, tmpY, '-', ...
        'Color', [1, 1, 1] .* sum(bValid(1:iC)) ./ (sum(bValid) + 1)); % different distances in different shades of gray
end
% plot COM
myColors    = circshift(hsv, size(hsv, 1) / 2);
potAngles   = linspace(-pi, pi, size(myColors, 1) + 1);
thisColor   = myColors(discretize(cfg.prefELRPD, potAngles), :);
plot(cfg.locdir_COMxy(1), cfg.locdir_COMxy(2), 'ko', ...
    'MarkerFaceColor', thisColor, 'LineWidth', 2);
% plot boundary
plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
t1 = text(0, 6000, 'x');
t2 = text(-6000, 0, 'y');
t3 = text(0, -6500, 'Reference point');
set([t1, t2, t3], ...
    'fontunits', 'centimeters', 'fontsize', 0.4, ...
    'HorizontalAlignment', 'center');
hold off;
set(gca, ...
    'ydir', 'reverse', ... % flip y-axis
    'xlim', [-5000, 5000], 'ylim', [-5000, 5000]);
axis off square equal;

%----- waveform
ax3 = axes('units', 'centimeters', 'position', [1, 4, 2, 2]);
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
xlabel('ms', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'units', 'normalized', 'position', [0.5, -0.1, 0]);
ylabel('\muV', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'units', 'normalized', 'position', [-0.12, 0.5, 0]);
text(0.5, 1.15, ['n=', num2str(cfg.nspk)], ...
    'units', 'normalized', ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, ...
    'HorizontalAlignment', 'center');
