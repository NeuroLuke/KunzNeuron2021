function out = LK_TH_locSmooth_170520(dt)
%
% LK_TH_locSmooth_170520 smoothes a locational firing rate map.
%
% Input is a structure with fields
% 	FR          --> firing-rate map
% 	smoothFac   --> smoothing factor
% 
% Output is a structure with fields
%   smFR        --> smoothed firing rate map
%
% Reference: Burgess et al., 2005
%
% Lukas Kunz, 2021

% unoccupied positions
bUnoccupied         = isnan(dt.FR);

% prepare smoothing
b                   = ones(dt.smoothFac, dt.smoothFac); % boxcar smooth function
c                   = ones(size(dt.FR));
c(bUnoccupied)    	= 0;
denom               = filter2(b, c);

% smooth firing rate map
FR                  = dt.FR;
FR(bUnoccupied)     = 0;
smFR                = filter2(b, FR);
smFR                = smFR ./ denom; % correct for amount of smoothing

% set unoccupied positions to nan
smFR(bUnoccupied)   = nan;

% create output
out         = [];
out.smFR    = smFR;

end