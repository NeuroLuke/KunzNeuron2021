function out = LK_TH_CompANOVA_210520(dt)
%
% LK_TH_CompANOVA_210520 computes an n-way ANOVA to estimate the influence
% of multiple predictors on neuronal firing rate.
% 
% Input is a structure with fields
%   FR                  --> firing rate per time bin (n x 1 vector)
%   group               --> level of each predictor per time bin (n x m matrix)
%   groupNames          --> predictor names (m x 1 cell)
%   ANOVA               --> ANOVA settings
%   bSkipA1 (optional)  --> whether to skip the computation of additional one-way ANOVAs
%
% Output is a structure with fields
%   groupNames          --> predictor names (m x 1 cell)
%   p                   --> p-values from the n-way ANOVA (m x 1 vector)
%   tbl                 --> statistics table from the n-way ANOVA
%   stats               --> statistics from the n-way ANOVA
%   terms               --> terms estimated in the n-way ANOVA
%   levelsOfInterest    --> predictor levels used in the n-way ANOVA
%   m                   --> estimated means from the n-way ANOVA
%   a1_m (optional)     --> estimated means from the predictor-wise one-way ANOVAs
%
% Lukas Kunz, 2021

%% preprocessing

% ensure that you have enough observations per group level
levelsOfInterest    = cell(size(dt.groupNames)); % e.g., size = 3 x 1; each cell entry contains the unique group-levels for each group
numObsMask          = nan(size(dt.group)); % timepoint-specific mask; e.g., size = 10,000 x 3
for iGroup = 1:size(dt.group, 2)
    
    % unique levels for this group
    bNan            = isnan(dt.group(:, iGroup));
    uniqueLevels    = unique(dt.group(~bNan, iGroup));
    
    % number of discrete observations per group level
    numObs          = nan(size(uniqueLevels, 1), 1);
    for iLevel = 1:size(uniqueLevels, 1)
        
        % get the discrete observations
        [~, NUM]            = bwlabel(dt.group(:, iGroup) == uniqueLevels(iLevel, 1));
        
        % count observations per level
        numObs(iLevel, 1)   = NUM;
    end
    
    % exclude levels for which the number of observations is too low
    levelsOfInterest{iGroup}    = uniqueLevels(numObs >= dt.ANOVA.minNumObsPerBin);
    numObsMask(:, iGroup)       = ismember(dt.group(:, iGroup), levelsOfInterest{iGroup});

    % cave: for specific predictors, do not apply a threshold regarding the
    % discrete number of observations
    if isfield(dt.ANOVA, 'excFromNumObs') && any(strcmp(dt.groupNames{iGroup}, dt.ANOVA.excFromNumObs)) % exceptions
        levelsOfInterest{iGroup}    = uniqueLevels;
        numObsMask(:, iGroup)       = ismember(dt.group(:, iGroup), levelsOfInterest{iGroup});
    end
end

% create mask for data: enough observations and no nans
bMask   = all(numObsMask, 2) & all(~isnan(dt.group), 2);

%% compute ANOVAN

% n-way ANOVA
[p, tbl, stats, terms]  = anovan(dt.FR(bMask), dt.group(bMask, :), ...
    'varnames', dt.groupNames, ...
    'display', 'off', ...
    'sstype', dt.ANOVA.sstype, ...
    'model', dt.ANOVA.model);

%% obtain adjusted means for each term via multcompare

% adjusted means for each factor
m   = cell(numel(p), 1);
for ip = 1:size(m, 1)
    [~, m{ip, 1}]   = multcompare(stats, 'display', 'off', 'dimension', ip);
end

%% create output

out                     = [];
out.groupNames          = dt.groupNames;
out.p                   = p;
out.tbl                 = tbl;
out.stats               = stats;
out.terms               = terms;
out.levelsOfInterest    = levelsOfInterest;
out.m                   = m; % estimated means

%% return here if desired

if isfield(dt, 'bSkipA1') && strcmp(dt.bSkipA1, 'yes')
    return;
end

%% compute one-way ANOVA separately for all factors to get uncorrected firing rate maps (not controlling for other factors)

% original means for each factor
a1_m    = cell(numel(p), 1);
for ip = 1:size(a1_m, 1)
    
    % one-way ANOVA
    [~, ~, a1_stats]    = anova1(dt.FR(bMask), dt.group(bMask, ip), 'off');
    
    % estimated means
    [~, a1_m{ip, 1}]    = multcompare(a1_stats, 'display', 'off');
end

%% add output

% original means
out.a1_m                = a1_m;

end