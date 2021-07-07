function out = LK_SmoothBearingDistanceMap_20201008(dt)
%
% LK_SmoothBearingDistanceMap_20201008 smoothes a bearing-distance
% firing-rate map.
%
% Cave: bearing is a circular variable, so the smoothing has to account for
% this.
%
% Reference: Burgess et al., 2005.
%
% Lukas Kunz, 2021

% sanity check: is the smoothing factor odd?
if mod(dt.smoothFac, 2) == 0
    error('The smoothing factor is not an odd number.')
end

% extent firing-rate map at both ends of bearing
extFac  = (dt.smoothFac - 1) / 2; % extension factor
extFR   = [dt.FR(:, end - extFac + 1:end), dt.FR, dt.FR(:, 1:extFac)]; % extended firing-rate map

% unoccupied positions
bUnoccupied     = isnan(extFR);

% smoothing kernel
b   = ones(dt.smoothFac, dt.smoothFac);
if isfield(dt, 'smoothType') && strcmp(dt.smoothType, 'gaussian')
    h  	= fspecial(dt.smoothType, dt.smoothFac, 2); % type, size, SD
    b   = b .* h;
end

% prepare smoothing
c                   = ones(size(extFR));
c(bUnoccupied)    	= 0;
denom               = filter2(b, c); % how much smoothing at each bin

% smooth firing rate map
FR                  = extFR;
FR(bUnoccupied)     = 0; % set to zero so that it does not affect neighboring bins
smFR                = filter2(b, FR);
smFR                = smFR ./ denom; % correct for amount of smoothing at each bin

% set unoccupied positions to nan
smFR(bUnoccupied)   = nan;

% crop to correct size
smFR                = smFR(:, extFac + 1:end - extFac);

% create output
out         = [];
out.smFR    = smFR;

end