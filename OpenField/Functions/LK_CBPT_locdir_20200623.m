function out = LK_CBPT_locdir_20200623(dt)
%
% LK_CBPT_locdir_20200623 evaluates the overall statistical significance of
% egocentric local reference point direction modulation maps.
% 
% Input is a structure with fields
%   locDir                  --> general settings
%   locRefPoints            --> local reference points
%   locRefPoints2Exclude    --> local reference points to exclude from the
%                               analysis
%   locdir_F                --> ANOVA F values for ELRPD tuning
%   locdir_F_surro          --> surrogate F values for ELRPD tuning
%   locdir_Frank            --> ANOVA F ranks for ELRPD tuning
%   locdir_corrFR           --> firing rates for ELRPD tuning
%   numSurrogates           --> number of surrogates
%
% Output is a structure with fields
%   maxSumF_emp             --> empirical statistic
%   maxSumF_surro           --> surrogate statistics
%   locdir_clusRank         --> rank of empirical statistic within
%                               surrogate statistics
%   COMxy                   --> reference point
%   COMxyClosestLRP         --> LRP closest to the reference point
%   largestCluster          --> reference field
%
% Lukas Kunz, 2021

%% significance threshold

sigThresh   = 0.95; % corresponding to alpha = 0.05

%% mask LRPs that shall be excluded from the analysis

% F-values
F                                           = dt.locdir_F;
F(dt.locRefPoints2Exclude == true)          = nan; % exclude non-candidate LRPs

% surrogate-F-values
F_surro                                     = dt.locdir_F_surro;
F_surro(dt.locRefPoints2Exclude == true, :) = nan;

% F-ranks
Frank                                       = dt.locdir_Frank;
Frank(dt.locRefPoints2Exclude == true)      = nan;

%% get data into 2D shape
% you order the data such that negative y-values and negative x-values are
% at the lower left corner of the matrix (positive y-values and positive
% x-values are at the upper right corner of the matrix)

% 2D maps of F-values, surrogate-F-values, and F-rank-values
FMap     	= flipud(reshape(F, numel(dt.locDir.yCenters), numel(dt.locDir.xCenters))); % rows refer to y-values; columns refer to x-values
FSurroMap   = flipud(reshape(F_surro, numel(dt.locDir.yCenters), numel(dt.locDir.xCenters), dt.numSurrogates));
FRankMap    = flipud(reshape(Frank, numel(dt.locDir.yCenters), numel(dt.locDir.xCenters)));

%% empiricals

% identify clusters of significant F-ranks and sum-up F-rank-values
% within these clusters
[L_emp, NUM_emp]    = bwlabel(FRankMap > sigThresh, 4); % note: the "4" indicates "conn", which can be 4 | 8 (default)
sumF                = nan(NUM_emp, 1); % summary statistic
for iNUM = 1:NUM_emp
    sumF(iNUM, 1)   = sum(FRankMap(L_emp == iNUM));
end
[maxSumF_emp, maxIdx_emp]   = max(sumF);

% if no maximum could be obtained
if isempty(maxSumF_emp)
    maxSumF_emp     = nan;
    largestCluster  = false(size(FMap));
else
    % location of largest cluster (= reference field)
    largestCluster  = L_emp == maxIdx_emp;
end

%% surrogates

% preallocate
surroIdx     	= 1:dt.numSurrogates;
maxSumF_surro 	= nan(dt.numSurrogates, 1);
for iSurro = 1:dt.numSurrogates
    
    % surrogate 2D-maps
    tmpFMap         = FSurroMap(:, :, surroIdx == iSurro); % one surrogate map
    tmpFSurroMap    = cat(3, FMap, FSurroMap(:, :, surroIdx ~= iSurro)); % empirical map and all other surrogate maps
    tmpFRankMap     = sum(tmpFMap > tmpFSurroMap, 3) ./ sum(~isnan(tmpFSurroMap), 3); % one surrogate rank map
    
    % identify clusters of significant F ranks and sum-up F-rank-values
    % within these clusters
    [L_surro, NUM_surro]    = bwlabel(tmpFRankMap > sigThresh, 4);
    sumF                    = nan(NUM_surro, 1); % summary statistic
    if NUM_surro > 0
        for iNUM = 1:NUM_surro
            sumF(iNUM, 1)           = sum(tmpFRankMap(L_surro == iNUM));
        end
        maxSumF_surro(iSurro, 1)    = max(sumF);
    end
end

% get cluster rank to assess overall significance
locdir_clusRank = sum(maxSumF_emp > maxSumF_surro) ./ sum(~isnan(maxSumF_surro));

%% identify center of mass of largest cluster (= reference point)

% center of mass of largest cluster, weighted by the F-ranks ("weCe")
stats       = regionprops(largestCluster, FRankMap, 'WeightedCentroid'); % regionprops(binary, weights, method)
if isempty(stats)
    weCe    = nan(1, 2);
else
    weCe    = stats.WeightedCentroid; % first element of WeightedCentroid is the horizontal coordinate (or x-coordinate)
end

% x-coordinate of COM
potIdx      = 1:numel(dt.locDir.xCenters);
xRatio      = (weCe(1) - min(potIdx)) / (max(potIdx) - min(potIdx));
COMx        = xRatio * range(dt.locDir.xCenters) - range(dt.locDir.xCenters) / 2;

% y-coordinate of COM (cave: y-values are flipped, thus subtract from 1)
potIdx      = 1:numel(dt.locDir.yCenters);
yRatio      = 1 - (weCe(2) - min(potIdx)) / (max(potIdx) - min(potIdx));
COMy        = yRatio * range(dt.locDir.yCenters) - range(dt.locDir.yCenters) / 2;

% get center of mass
COMxy       = [COMx, COMy];

%% closest xy-bin of the largest cluster towards the center of mass

% unfold largest cluster, so that it matches the sequence of LRPs
unfLargestCluster       = reshape(flipud(largestCluster), numel(dt.locDir.xCenters) * numel(dt.locDir.yCenters), 1);

% xy coordinates of LRPs inside the largest cluster
LRPwithinLargestCluster = dt.locRefPoints(unfLargestCluster, :);

% calculate distance of center of mass to all LRPs within the largest
% cluster
D                       = pdist2(COMxy, LRPwithinLargestCluster);
[~, minDidx]            = min(D);

% closest LRP
COMxyClosestLRP     	= LRPwithinLargestCluster(minDidx, :);

%% create output

out                   	= [];
% cluster summary statistic and rank
out.maxSumF_emp        	= maxSumF_emp;
out.maxSumF_surro      	= maxSumF_surro;
out.locdir_clusRank    	= locdir_clusRank;
% center of mass
out.COMxy              	= COMxy;
out.COMxyClosestLRP   	= COMxyClosestLRP;
% unfolded largest cluster
out.largestCluster   	= unfLargestCluster;

