function yl = myyline(y, myLineStyle, myColor, myStackLevel)
%
% MYYLINE creates a constant line at a specified x value.
%
% MYYLINE is similar to MATLAB's yline, but allows to stack the line to the
% bottom.
%
% Lukas Kunz, 2021

% information about current axes
tmpAx = get(gca);

% plot yline
yl = plot([min(tmpAx.XLim), max(tmpAx.XLim)], [y, y], myLineStyle, ...
    'Color', myColor);

% stack
uistack(yl, myStackLevel);