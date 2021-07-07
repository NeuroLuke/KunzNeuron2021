function out = LK_EvalBehOverall_20200623(dt)
%
% LK_EvalBehOverall_20200623 processes the behavioral data of the
% OpenField task.
%
% Input is a structure with fields
%   paths               --> path information
%   subject             --> subject information
%   direc (optional)    --> structure with settings for extracting direction information
%   place (optional)    --> structure with settings for extracting place information
%   locDir (optional)   --> structure with settings for extracting egocentric local direction information
%   locDist (optional)  --> structure with settings for extracting egocentric local distance information
%   maxAbsYaw           --> maximum absolute yaw value (based on the original virtual arena settings; default: 32768)
%   maxNumTrials        --> maximum number of trials to include (default: 'all')
%   newTimeRes          --> temporal resolution for new timeline (default: 0.1)
%   naviCutoff          --> cutoff for navigation periods (default: 0.001)
%   naviSmooth          --> setting for smoothing the navigation-periods vector (default: 2)
%   rotSmooth           --> setting for smoothing the rotation-periods vector (default: 2)
%   ...
%
% Output is a structure with fields
%   orig                --> various behavioral information along original timeline
%   new                 --> various behavioral information along new timeline
%
% Lukas Kunz, 2021

fprintf('\nBehavioral processing of %s ...\n', dt.subject);

%% load and adjust data

% load "trials", which contains information about trialidx, objectidx, ITI,
% cue, retrieval, feedback, reencoding, grab, correctXY, dropXY, drop error
tmp     = load(strcat(dt.paths.beh, dt.subject, filesep, 'trials.mat'));
trials  = tmp.trials;
fprintf('Original size of "trials": \t%d x %d.\n', size(trials));

% restrict to a specific number of trials
if ~strcmp(dt.maxNumTrials, 'all')
    
    % restrict to complete trials
    if strcmp(dt.maxNumTrials, 'completeTrialsOnly')
        if size(trials, 1) > 160        
            trials  = trials(1:160, :); % restrict to 160 complete trials
        else
            trials  = trials(1:size(trials, 1) - 1, :); % remove incomplete (i.e., last) trial
        end
        
    % restrict to specific number of trials    
    elseif size(trials, 1) > dt.maxNumTrials
        trials  = trials(1:dt.maxNumTrials, :);
    end
end
fprintf('New size of "trials": \t%d x %d.\n', size(trials));

% load "behinfo" that contains, for each timepoint, information about time,
% xyz, yaw, orientation, trialidx, trialphase, objectidx, correctXY,
% dropXY, drop error
tmp     = load(strcat(dt.paths.beh, dt.subject, filesep, 'behinfo.mat'));
behinfo = tmp.behinfo;
fprintf('Original size of "behinfo": \t%d x %d.\n', size(behinfo));

% restrict analysis period to period between first and last cue (this
% is also the time period during which spikes were extracted)
twoi        = [trials(1, 4), trials(end, 4)];
if strcmp(dt.maxNumTrials, 'completeTrialsOnly')
    warning('Working with complete trials only.');
    twoi    = [trials(1, 4), trials(end, 8)];
end
if isnan(trials(end, 4)) % sanity check
    error('Last cue timepoint is a NaN.');
end
behinfo     = behinfo(behinfo(:, 1) >= twoi(1) & behinfo(:, 1) <= twoi(2), :); % cut behavioral information to time window of interest
fprintf('New size of "behinfo": \t%d x %d.\n', size(behinfo));
fprintf('The data segment lasts %.2f min.\n', range(behinfo(:, 1)) / 60);

% object index per trial
objIdx  = trials(:, 2);

%% ======================================================================== ORIGINAL TIMELINE

%% get behavioral information

% store entire behavioral information
orig                = [];
orig.behinfo        = behinfo;

% estimate further information
orig.time           = orig.behinfo(:, 1);
orig.durations      = [diff(orig.time); median(diff(orig.time))]; % add a median duration at the end
orig.distances      = [sqrt(diff(orig.behinfo(:, 2)) .^ 2 + diff(orig.behinfo(:, 3)) .^ 2); 0]; % add a zero-distance at the end (no moving)
orig.speed          = orig.distances ./ orig.durations;
orig.speedchange    = [0; diff(orig.speed)]; % add a zero-speed-change at the beginning
orig.acceleration   = orig.speedchange ./ orig.durations;
orig.trialPhase     = orig.behinfo(:, 8); % trial phase: 1 = ITI, 2 = cue, 3 = retrieval, 4 = feedback, 5 = reencoding
orig.bActive        = orig.trialPhase > 2; % subject is active during retrieval/feedback/re-encoding

% navigation periods
orig.bNavi          = orig.speed > dt.naviCutoff;
% smooth navigation periods
if isfield(dt, 'naviSmooth') && dt.naviSmooth > 0
    convVec      	= ones(dt.naviSmooth / 0.02, 1); % convolution vector; original temporal resolution is ~20 msec
    orig.bNavi   	= conv(orig.bNavi, convVec, 'same') > 0;
end
% sanity constraint: restrict navigation periods to
% retrieval/feedback/reencoding
orig.bNavi          = orig.bNavi & orig.bActive;

%% direction

if isfield(dt, 'direc')
    
    % report
    fprintf('Processing direction information.\n');
    
    % get yaws and discretize yaws into yaw bins
    orig.yaws       = (orig.behinfo(:, 5) ./ dt.maxAbsYaw) .* pi; % map original yaw-values onto [-pi, pi]
    orig.yawBin     = discretize(orig.yaws, dt.direc.angleEdges);
    
    % calculate navigation direction based on subsequent xy locations
    % cave: navigation direction is similar, but not identical, to yaws
    % and can only be determined for navigation periods
    orig.navDir     = [atan2(diff(orig.behinfo(:, 3)), diff(orig.behinfo(:, 2))); 0]; % atan2(y, x); add a 0° navigation direction at the end
    orig.navDirBin  = discretize(orig.navDir, dt.direc.angleEdges);
    
    % estimate angular speed and acceleration
    orig.angDistances       = abs([angdiff(orig.yaws); 0]); % add a zero-distance at the end (no turning)
    orig.angSpeed       	= orig.angDistances ./ orig.durations;
    orig.angSpeedChange     = [0; diff(orig.angSpeed)]; % add a zero-angular-speed-change at the beginning
    orig.angSpeedChange(abs(orig.angSpeedChange) > quantile(abs(orig.angSpeedChange), 0.95))    = 0; % remove angular speed changes that are too fast (implausible)
    orig.angAcceleration    = orig.angSpeedChange ./ orig.durations;
    
    % rotation mask
    orig.bRotation          = orig.angSpeed > 0;
    % smooth rotation mask
    if isfield(dt, 'rotSmooth') && dt.rotSmooth > 0
        convVec             = ones(dt.rotSmooth / 0.02, 1); % convolution vector; original temporal resolution is ~20 msec
        orig.bRotation      = conv(orig.bRotation, convVec, 'same') > 0;
    end
    % restrict to periods of retrieval/feedback/reencoding
    orig.bRotation          = orig.bRotation > 0 & orig.bActive;
end

%% place

if isfield(dt, 'place')
    
    % report
    fprintf('Processing place information.\n');
    
    % discretize locations into locational bins
    orig.xBin 	= discretize(orig.behinfo(:, 2), dt.place.xEdges);
    orig.yBin 	= discretize(orig.behinfo(:, 3), dt.place.yEdges);
    orig.xyBin 	= (orig.xBin - 1) .* dt.place.locRes + orig.yBin; % combined xy-bin
end

%% egocentric direction towards local reference point

if isfield(dt, 'locDir')
    
    % report
    fprintf('Processing egocentric direction to local reference point [x = %.1f, y = %.1f].\n', ...
        dt.locDir.refPoint(1), dt.locDir.refPoint(2));
    
    % compute direction between subject position and local reference point
    orig.alloLocDir     = atan2(dt.locDir.refPoint(2) - orig.behinfo(:, 3), dt.locDir.refPoint(1) - orig.behinfo(:, 2)); % atan2(y, x)
    
    % compute egocentric direction to local reference point
    orig.egoLocDir      = angdiff(orig.yaws, orig.alloLocDir); % angdiff(x, y) returns the angular difference delta = y - x
    orig.egoLocDirBin   = discretize(orig.egoLocDir, dt.locDir.angleEdges);
end

%% egocentric distance towards local reference point

if isfield(dt, 'locDist')

    % report
    fprintf('Processing egocentric distance to local reference point [x = %.1f, y = %.1f].\n', ...
        dt.locDist.refPoint(1), dt.locDist.refPoint(2));
    
    % compute distance between the subject's current location and the
    % reference location
    orig.egoLocDist     = transpose(pdist2(dt.locDist.refPoint, orig.behinfo(:, 2:3)));
    orig.egoLocDistBin  = discretize(orig.egoLocDist, dt.locDist.bin_edges);
end

%% ======================================================================== NEW TIMELINE

fprintf('Creating a new timeline with a temporal resolution of %.1f Hz.\n', 1 / dt.newTimeRes);

%% create new timeline

% define new timeline
new                 = [];
new.timeRes         = dt.newTimeRes;
newTimeLine         = transpose(min(orig.time):new.timeRes:max(orig.time));

% create new "behinfo"-file
new.behinfo         = nan(size(newTimeLine, 1), size(orig.behinfo, 2));
new.behinfo(:, 1)   = newTimeLine;

% preallocate other variables
if isfield(dt, 'direc')
    new.yaws            = nan(size(newTimeLine, 1), 1);
    new.yawBin          = nan(size(newTimeLine, 1), 1);
    new.navDir          = nan(size(newTimeLine, 1), 1);
    new.navDirBin       = nan(size(newTimeLine, 1), 1);
end
if isfield(dt, 'locDir')
    new.egoLocDir       = nan(size(newTimeLine, 1), 1);
    new.egoLocDirBin    = nan(size(newTimeLine, 1), 1);
end
if isfield(dt, 'locDist')
    new.egoLocDist      = nan(size(newTimeLine, 1), 1);
    new.egoLocDistBin   = nan(size(newTimeLine, 1), 1);
end

% behinfo
for iT = 1:size(new.behinfo, 1) - 1
    
    % get section of original behinfo-file with data from new time window
    bTWOI           = orig.behinfo(:, 1) >= new.behinfo(iT, 1) & orig.behinfo(:, 1) < new.behinfo(iT + 1, 1);
    thisBehinfo     = orig.behinfo(bTWOI, :); % time, xyz, yaw, orientation, trialidx, trialphase, objectidx, correctXY, dropXY, drop error
    
    %% behinfo
    
    % x, y, z
    new.behinfo(iT, 2:4)            = mean(thisBehinfo(:, 2:4), 1); % average along first dimension
    
    % yaw, orientation (angular variables)
    thisYaws                        = (thisBehinfo(:, 5) ./ dt.maxAbsYaw) .* pi; % convert onto [-pi, pi] so that you can calculate the circular mean
    new.behinfo(iT, 5:6)            = [(circ_mean(thisYaws) ./ pi) .* dt.maxAbsYaw, nan]; % map back into original "yaw-space"; don't fill in a value for orientation, because not necessary  
    
    % trialidx, trialphase, objectidx, correctXY, dropXY, error
    new.behinfo(iT, 7:14)         	= mode(thisBehinfo(:, 7:14), 1); % along first dimension
    
    %% direction
    
    if isfield(dt, 'direc')
        
        % yaws
        new.yaws(iT, 1)         	= circ_mean(orig.yaws(bTWOI)); % already within [-pi, pi]
        new.yawBin(iT, 1)        	= mode(orig.yawBin(bTWOI), 1);
        
        % navigation direction
        new.navDir(iT, 1)       	= circ_mean(orig.navDir(bTWOI));
        new.navDirBin(iT, 1)      	= mode(orig.navDirBin(bTWOI), 1);
    end
    
    %% egocentric direction to local reference point, new timeline
    
    if isfield(dt, 'locDir')
        
        % egocentric direction
        new.egoLocDir(iT, 1)     	= circ_mean(orig.egoLocDir(bTWOI));
        new.egoLocDirBin(iT, 1)    	= mode(orig.egoLocDirBin(bTWOI), 1);
    end
    
    %% egocentric distance to local reference point, new timeline
    
    if isfield(dt, 'locDist')
        
        % egocentric distance
        new.egoLocDist(iT, 1)       = mean(orig.egoLocDist(bTWOI));
        new.egoLocDistBin(iT, 1)    = mode(orig.egoLocDistBin(bTWOI), 1);
    end        
end

% additional behavioral information
new.time            = new.behinfo(:, 1);
new.durations       = [diff(new.time); median(diff(new.time))]; % add a median duration at the end
new.distances       = [sqrt(diff(new.behinfo(:, 2)) .^ 2 + diff(new.behinfo(:, 3)) .^ 2); 0]; % add a zero-distance at the end (no moving)
new.speed           = new.distances ./ new.durations;
new.speedchange     = [0; diff(new.speed)]; % add a zero-speed-change at the beginning
new.acceleration    = new.speedchange ./ new.durations;
new.trialPhase      = new.behinfo(:, 8);
new.bActive         = new.trialPhase > 2; % subject is active during retrieval/feedback/re-encoding

% navigation mask
new.bNavi           = new.speed > dt.naviCutoff;
% smooth navigation mask
if isfield(dt, 'naviSmooth') && dt.naviSmooth > 0
    convVec     	= ones(dt.naviSmooth / dt.newTimeRes, 1); % convolution vector
    new.bNavi      	= conv(new.bNavi, convVec, 'same') > 0;
end
% sanity constraint: restrict navigation periods to retrieval/feedback/reencoding
new.bNavi           = new.bNavi & new.bActive;

% direction
if isfield(dt, 'direc')
    
    % estimate angular speed and acceleration
    new.angDistances        = abs([angdiff(new.yaws); 0]); % add a zero-distance at the end (no turning)
    new.angSpeed            = new.angDistances ./ new.durations;
    new.angSpeedChange      = [0; diff(new.angSpeed)]; % add a zero-angular-speed-change at the beginning
    new.angSpeedChange(abs(new.angSpeedChange) > quantile(abs(new.angSpeedChange), 0.95))   = 0; % remove angular speed changes that are too fast (implausible)
    new.angAcceleration     = new.angSpeedChange ./ new.durations;
    
    % rotation mask
    new.bRotation           = new.angSpeed > 0;
    % smooth rotation mask
    if isfield(dt, 'rotSmooth') && dt.rotSmooth > 0
        convVec             = ones(dt.rotSmooth / dt.newTimeRes, 1); % convolution vector
        new.bRotation       = conv(new.bRotation, convVec, 'same') > 0;
    end
    % restrict to periods of retrieval/feedback/reencoding
    new.bRotation           = new.bRotation > 0 & new.bActive;
end

% place
if isfield(dt, 'place')
    
    % discretize locations into locational bins
    new.xBin  	= discretize(new.behinfo(:, 2), dt.place.xEdges);
    new.yBin    = discretize(new.behinfo(:, 3), dt.place.yEdges);
    new.xyBin   = (new.xBin - 1) .* dt.place.locRes + new.yBin;
end

%% collect output

% initialize
out                 = [];

% data from original timeline
out.orig            = orig;
out.orig.trials   	= trials;
out.orig.objIdx    	= objIdx;

% data from new timeline
out.new             = new;
out.new.trials   	= trials;
out.new.objIdx     	= objIdx;

