function out = LK_EstimateBearingDistanceMap_20210515(cfg)
%
% LK_EstimateBearingDistanceMap_20210515 estimates egocentric 2D
% bearing-distance maps.
%
% Input is a structure with fields
%   FR                  --> firing rates
%   dist2COM            --> distances to the reference point
%   egoLocDir           --> egocentric directions to the reference point
%   bearingBinEdges    	--> edges for bearing bins
%   bearingBinCenters   --> centers for bearing bins
%   distanceBinEdges  	--> edges for distance bins
%   distanceBinCenters  --> centers for distance bins
%   bSmooth             --> boolean whether to smooth the map
%   smoothFac           --> smoothing factor
%   smoothType          --> type of kernel for smoothing
%
% Output is a structure with fields
%   firingRate          --> bearing-distance map of (smoothed) firing rates
%   meanBearingTuning   --> average bearing tuning
%   meanDistanceTuning  --> average distance tuning
%
% Lukas Kunz, 2021

%% firing rate as a function of distance and bearing

% bin behavioral data
distBin         = discretize(cfg.dist2COM, cfg.distanceBinEdges);
egoLocDirBin  	= discretize(cfg.egoLocDir, cfg.bearingBinEdges);

% create firing-rate map
firingRate   	= nan(numel(cfg.distanceBinCenters), numel(cfg.bearingBinCenters)); % y: distance; x: bearing
for iB = 1:numel(cfg.bearingBinCenters)
    for iD = 1:numel(cfg.distanceBinCenters)
        firingRate(iD, iB) 	= mean(cfg.FR(egoLocDirBin == iB & distBin == iD)); % top to bottom: small to large distances
    end
end

% smooth
if isfield(cfg, 'bSmooth') && cfg.bSmooth == true
    
    % smooth using a 2D kernel
    dt              = [];
    dt.FR           = firingRate;
    dt.smoothFac    = cfg.smoothFac;
    dt.smoothType   = cfg.smoothType;
    outSmooth       = LK_SmoothBearingDistanceMap_20201008(dt);
    
    % replace by smoothed firing-rate map
    firingRate   	= outSmooth.smFR;
end

%% mean bearing tuning and mean distance tuning

% means across bins
meanBearingTuning   = nanmean(firingRate, 1);
meanDistanceTuning  = nanmean(firingRate, 2);

%% collect output

out                     = [];
out.firingRate          = firingRate;
out.meanBearingTuning  	= meanBearingTuning;
out.meanDistanceTuning  = meanDistanceTuning;

