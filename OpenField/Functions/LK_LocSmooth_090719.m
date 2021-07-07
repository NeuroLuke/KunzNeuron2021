function out = LK_LocSmooth_090719(dt)
%
% LK_LocSmooth_090719 smoothes a locational firing rate map.
% The smoothing factor (smoothFac) should be an odd number.
%
% Reference: Burgess et al., 2005
%
% Lukas Kunz, 2021

% unoccupied positions
bUnoccupied         = isnan(dt.FR);

% prepare smoothing
b                   = ones(dt.smoothFac, dt.smoothFac);
c                   = ones(size(dt.FR));
c(bUnoccupied)    	= 0;
denom               = filter2(b, c);

% smooth firing rate map
FR                  = dt.FR;
FR(bUnoccupied)     = 0; % set to zero so that it does not affect neighboring bins
smFR                = filter2(b, FR);
smFR                = smFR ./ denom; % correct for amount of smoothing

% set unoccupied positions to nan
smFR(bUnoccupied)   = nan;

% create output
out         = [];
out.smFR    = smFR;

end