function xl = myxline(x, myLineStyle, myColor, myStackLevel)
%
% MYXLINE creates a constant line at a specified x value.
%
% MYXLINE is similar to MATLAB's xline, but allows to stack the line to the
% bottom.
%
% Lukas Kunz, 2021

% information about current axes
tmpAx = get(gca);

% plot xline
xl = plot([x, x], [min(tmpAx.YLim), max(tmpAx.YLim)], myLineStyle, ...
    'Color', myColor);

% stack
uistack(xl, myStackLevel);