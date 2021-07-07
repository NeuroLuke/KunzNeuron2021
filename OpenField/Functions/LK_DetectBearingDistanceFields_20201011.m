function out = LK_DetectBearingDistanceFields_20201011(cfg)
%
% LK_DetectBearingDistanceFields_20201011 detects firing fields in 2D
% bearing-distance firing-rate maps.
%
% Input: structure with fields
%   - firingRate: empirical 2D firing-rate map
%   - firingRateSurro: n surrogate 2D firing-rate maps
%
% Lukas Kunz, 2021

% detect fields with increased firing rates
firingRateRanks = sum(cfg.firingRate > cfg.firingRateSurro, 3) ./ sum(~isnan(cfg.firingRateSurro), 3);
bField          = firingRateRanks > 0.95; % uncorrected alpha-level of 5%
[L, NUM]        = bwlabel(bField, 4); % "4" denotes the connectivity

% correct for circularity of the bearing angles (using a connectivity of 4)
for iDistance = 1:size(L, 1)
    if bField(iDistance, end) == true && bField(iDistance, 1) == true
        L(L == L(iDistance, end))  = L(iDistance, 1); % "relabel"
    end
end

% estimate field size or field strength
fieldSize   = nan(NUM, 1);
for iNUM = 1:NUM
    fieldSize(iNUM, 1)  = sum(sum(L == iNUM));
end

% identify largest field and use it as the bearing-distance field
[maxFieldSize, idxMaxFieldSize]     = max(fieldSize);

% if there is a field or if there is no field
if ~isempty(maxFieldSize)
    bearingDistanceField    = L == idxMaxFieldSize;
elseif isempty(maxFieldSize)
    fieldSize               = nan;
    maxFieldSize            = nan;
    bearingDistanceField    = bField;
end

% collect output
out                         = [];
out.bField                  = bField;
out.fieldSize               = fieldSize;
out.maxFieldSize            = maxFieldSize;
out.bearingDistanceField    = bearingDistanceField;
