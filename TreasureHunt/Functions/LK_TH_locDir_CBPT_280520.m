function out = LK_TH_locDir_CBPT_280520(dt)
%
% LK_TH_locDir_CBPT_280520 evaluates the overall statistical significance
% of egocentric local reference point direction modulation maps.
%
% Input is a structure with fields
%   locDir              --> general settings
%   all_locdir_F        --> ANOVA F values for ELRPD tuning
%   all_locdir_Fsurro   --> surrogate F values for ELRPD tuning
%   all_locdir_Frank    --> ANOVA F ranks for ELRPD tuning
%   numSurrogates       --> number of surrogates
%
% Output is a structure with fields
%   maxSumStat_emp      --> empirical statistic
%   locdir_clusRank     --> rank of empirical within surrogate statistics
%   COMxz               --> reference point
%   COMxzClosestLRP     --> LRP closest to the reference point
%   largestCluster      --> reference field
% 
% Lukas Kunz, 2021

%% significance threshold

sigThresh   = 0.95; % corresponding to alpha = 0.05

%% mask LRPs that shall be excluded from the analysis

% F-values
F                                           = dt.all_locdir_F;
F(dt.locDir.LRPs2Exclude == true)           = nan; % exclude non-candidate LRPs

% surrogate-F-values
Fsurro                                      = dt.all_locdir_Fsurro;
Fsurro(dt.locDir.LRPs2Exclude == true, :)   = nan;

% F-ranks
Frank                                       = dt.all_locdir_Frank;
Frank(dt.locDir.LRPs2Exclude == true)       = nan;

%% get data into 2D shape
% you order the data such that high z-values and high x-values are at the
% lower left corner of the matrix (low z-values and low x-values are at the
% upper right corner of the matrix)

% 2D maps of F-values, surrogate-F-values, and F-rank-values
FMap     	= flipud(reshape(F, numel(dt.locDir.zCenters), numel(dt.locDir.xCenters)));
FSurroMap   = flipud(reshape(Fsurro, numel(dt.locDir.zCenters), numel(dt.locDir.xCenters), dt.numSurrogates));
FRankMap    = flipud(reshape(Frank, numel(dt.locDir.zCenters), numel(dt.locDir.xCenters)));

%% empiricals

% identify clusters of significant F-ranks and sum up the F-rank-values
% within these clusters
[L_emp, NUM_emp]    = bwlabel(FRankMap > sigThresh, 4); % note: the "4" indicates "conn", which can have 4 | 8 (default)
sumStat             = nan(NUM_emp, 1); % summary statistic
for iNUM = 1:NUM_emp
    sumStat(iNUM, 1)   = sum(FRankMap(L_emp == iNUM));
end
[maxSumStat_emp, maxIdx_emp]    = max(sumStat); % maximum summary statistic

% if no maximum is obtained
if isempty(maxSumStat_emp)
    maxSumStat_emp  = nan;
    largestCluster  = false(size(FMap));
else
    % store information about location of largest cluster
    largestCluster  = L_emp == maxIdx_emp;
end

%% surrogates

% preallocate
surroIdx            = 1:dt.numSurrogates;
maxSumStat_surro    = nan(dt.numSurrogates, 1);
for iSurro = 1:dt.numSurrogates
    
    % surrogate 2D maps
    tmpFMap         = FSurroMap(:, :, surroIdx == iSurro); % one surrogate map
    tmpFSurroMap    = cat(3, FMap, FSurroMap(:, :, surroIdx ~= iSurro)); % empirical map and all other surrogate maps
    tmpFRankMap     = sum(tmpFMap > tmpFSurroMap, 3) ./ sum(~isnan(tmpFSurroMap), 3); % one surrogate rank map
    
    % identify clusters of significant F ranks and sum-up F-rank-values
    % within these clusters
    [L_surro, NUM_surro]    = bwlabel(tmpFRankMap > sigThresh, 4); % note: the "4" indicates "conn", which can have 4 | 8 (default)
    sumStat               	= nan(NUM_surro, 1); % summary statistic
    if NUM_surro > 0
        for iNUM = 1:NUM_surro
            sumStat(iNUM, 1)        = sum(tmpFRankMap(L_surro == iNUM));
        end
        maxSumStat_surro(iSurro, 1) = max(sumStat);
    end
end

% cluster rank to assess overall significance
locdir_clusRank = sum(maxSumStat_emp > maxSumStat_surro) ./ sum(~isnan(maxSumStat_surro));

%% identify center of mass for largest cluster

% center of mass of largest cluster, weighted by the F ranks ("weCe")
stats      	= regionprops(largestCluster, FRankMap, 'WeightedCentroid'); % regionprops(binary, weights, method)
if isempty(stats)
    weCe  	= nan(1, 2);
else
    weCe  	= stats.WeightedCentroid;
end

% x-coordinate of COM
potIdx    	= 1:numel(dt.locDir.xCenters);
xRatio   	= (weCe(1) - min(potIdx)) / (max(potIdx) - min(potIdx));
COMx      	= max(dt.locDir.xCenters) - xRatio * range(dt.locDir.xCenters);

% z-coordinate of COM
potIdx     	= 1:numel(dt.locDir.zCenters);
zRatio    	= 1 - (weCe(2) - min(potIdx)) / (max(potIdx) - min(potIdx));
COMz       	= max(dt.locDir.zCenters) - zRatio * range(dt.locDir.zCenters);

% get center of mass
COMxz       = [COMx, COMz];

%% get closest xz-bin of the largest cluster towards the center of mass

% unfold largest cluster, so that it matches the sequence of LRPs
unfLargestCluster       = reshape(flipud(largestCluster), numel(dt.locDir.zCenters) * numel(dt.locDir.xCenters), 1);

% xz coordinates of LRPs inside the largest cluster
LRPwithinLargestCluster = dt.locDir.LRPs(unfLargestCluster, :);

% calculate distance of center of mass to all LRPs within the largest
% cluster
D                       = pdist2(COMxz, LRPwithinLargestCluster);
[~, minDidx]            = min(D);

% closest LRP
COMxzClosestLRP     	= LRPwithinLargestCluster(minDidx, :);

%% create output

out                 	= [];
out.maxSumStat_emp    	= maxSumStat_emp; % cluster summary statistic
out.locdir_clusRank   	= locdir_clusRank; % rank of cluster summary statistic
out.COMxz            	= COMxz; % center of mass
out.COMxzClosestLRP 	= COMxzClosestLRP; % closest LRP to center of mass
out.largestCluster    	= unfLargestCluster; % unfolded largest cluster
