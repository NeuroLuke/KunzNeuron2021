function [f, percPerReg] = LK_PlotCellsPerRegion(dt)
%
% LK_PlotCellsPerRegion creates a barplot to present the percentage of
% cells per region.
%
% Input is a structure with fields
%   figEx           --> extension of figure (cm)
%   uniqueRegions   --> unique regions under examination
%   allRegions      --> region of each cell
%   bCell           --> boolean whether cell belongs to specific cell type
%   ylabel          --> label of y-axis
%   
% Output:
%   f               --> figure handle
%   percPerReg      --> percentage of specific cell type per region
%
% Lukas Kunz, 2021

% create figure
f = figure('units', 'centimeters', 'position', dt.figEx);
hold on;

% loop through regions
percPerReg  = nan(numel(dt.uniqueRegions), 1); % percentage of specific cell type per region
for iReg = 1:numel(dt.uniqueRegions)
    
    % calculate number of neurons for this region
    bThisReg    = strcmp(dt.allRegions, dt.uniqueRegions(iReg));
    myNum       = sum(bThisReg);
    
    % calculate percentage and p-value for this region
    myPerc                  = 100 * sum(bThisReg & dt.bCell) / sum(bThisReg);
    percPerReg(iReg, 1)     = myPerc;
    myAlpha                 = 0.05; % alpha level
    myChancePerc            = 100 * myAlpha; % chance percentage
    myPVal                  = myBinomTest(sum(bThisReg & dt.bCell), sum(bThisReg), myAlpha);
    
    % plot bar
    b1 = bar(iReg, myPerc);
    set(b1, ...
        'FaceColor', [0.7 0.7 0.7]);
    % number of cells in this region
    text(b1.XData, 0, num2str(myNum), ...
        'FontUnits', 'centimeters', 'fontsize', 0.4, ...
        'Color', [1, 1, 1], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
    % indicate significance
    if myPerc > myChancePerc
        if myPVal < 0.001
            text(b1.XData, b1.YData, '***', ...
                'fontunits', 'centimeters', 'fontsize', 0.5, ...
                'verticalalignment', 'baseline', 'horizontalalignment', 'center')
        elseif myPVal < 0.01
            text(b1.XData, b1.YData, '**', ...
                'fontunits', 'centimeters', 'fontsize', 0.5, ...
                'verticalalignment', 'baseline', 'horizontalalignment', 'center')
        elseif myPVal < 0.05
            text(b1.XData, b1.YData, '*', ...
                'fontunits', 'centimeters', 'fontsize', 0.5, ...
                'verticalalignment', 'baseline', 'horizontalalignment', 'center')
        end
    end
end
hold off;

% indicate chance level
yline(myChancePerc, '--', ...
    'Color', [0, 0, 0], 'LineWidth', 1);

% enhance axes
set(gca, ...
    'xlim', [0.2, numel(dt.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(dt.uniqueRegions), 'xticklabel', dt.uniqueRegions, 'tickdir', 'out', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
xlabel('Region', ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
ylabel(dt.ylabel, ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
