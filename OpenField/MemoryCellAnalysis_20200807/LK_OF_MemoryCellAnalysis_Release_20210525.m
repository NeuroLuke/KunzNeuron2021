%==========================================================================
% This script examines memory cells and whether they overlap with ELRPD
% cells.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.myRNG             = 3333;
param.clusterName       = 'HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768; % maximum absolute yaw value
param.newTimeRes        = 0.1; % time resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates     = 101; % number of surrogates
param.trialPhase        = 'Active'; % trial phase to use
%- for direction
direc                   = [];
direc.angularRes        = 30; % angular resolution
direc.angleEdges        = deg2rad(-180:direc.angularRes:180);
direc.angleCenters      = transpose(movmean(direc.angleEdges, 2, 'endpoints', 'discard'));
%- for location
place                   = [];
place.locRes            = 10; % spatial resolution
place.idxTemplate       = flipud(reshape(1:place.locRes^2, place.locRes, place.locRes));
place.xEdges            = linspace(-4500, 4500, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.yEdges            = linspace(-4500, 4500, place.locRes + 1);
place.yCenters          = movmean(place.yEdges, 2, 'endpoints', 'discard');
%- for memory
memor                   = [];
memor.corrType          = 'Pearson'; % correlation type
memor.bRankData         = 'yes'; % whether to rank data before correlation
memor.tailAlpha         = 0.025; % alpha-level at both ends of the surrogate distributions
memor.tailThresh        = 1 - memor.tailAlpha; % rank-threshold at both ends of the surrogate distributions

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectInformation_220318\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh       = 'E:\OpenField\Beh_210318\';
paths.save      = strcat('E:\OpenField\MemoryCellAnalysis_20200807\20201128_', param.trialPhase, '_', ...
    num2str(param.myRNG), '_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
    '_Memor', memor.corrType, '\', ...
    'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotation_Smooth', num2str(param.naviSmooth), '\');
if strcmp(memor.bRankData, 'yes')
    paths.save  = strcat('E:\OpenField\MemoryCellAnalysis_20200807\20201128_', param.trialPhase, '_', ...
        num2str(param.myRNG), '_Dir', num2str(direc.angularRes), ...
        '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
        '_RankedMemor', memor.corrType, '\', ...
        'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotation_Smooth', num2str(param.naviSmooth), '\');
end
mkdir(paths.save);
paths.arena     = 'E:\OpenField\Arena\';

% subjects
subjects    = {...
    'OF_001a'; ...
    'OF_001b'; ...
    'OF_002'; ...
    'OF_003a'; ...
    'OF_003b'; ...
    'OF_004a'; ...
    'OF_004b'; ...
    'OF_005'; ...
    'OF_006'; ...
    'OF_007'; ...
    'OF_008'; ...
    'OF_009'; ...
    'OF_010'; ...
    'OF_011'; ...
    'OF_012'; ...
    'OF_013a'; ...
    'OF_013b'; ...
    'OF_014'; ...
    };

%% save settings
save(strcat(paths.save, 'settings'));

%% random locations in the arena to convert drop errors into memory performances

% reset rng
rng(param.myRNG);

% many random locations within arena
dt              = [];
dt.maxR         = 5000; % maximum radius
dt.minR         = 0; % minimum radius
dt.N            = 5000001; % number of locations to create
dt.centerX      = 0;
dt.centerY      = 0;
locsInArena     = LK_RandomPointsInCircle_101119(dt);

%% preallocations

% load reference analysis with results on ELRPD cells
ELRPD   = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'locDir', 'allBeh');

% main results for each cell
allRes   	= [];

% for bookkeeping
exByWC      = []; % excluded by wave-clus
exByRefAna  = []; % excluded by reference analysis

% for behavioral control analyses
allBeh      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(param.myRNG);
    
    % subject-specific save path
    subjSavePath    = strcat(paths.save, subjects{iSub}, '\');
    if ~exist(subjSavePath, 'dir')
        mkdir(subjSavePath);
    end
    
    %% get interesting behavioral data
    
    % extract behavioral data
    behFile   	= strcat(subjSavePath, 'outBeh.mat');
    if exist(behFile, 'file') == 2
        
        % load previously saved behavioral data
        fprintf('Loading behavioral data from "%s".\n', behFile);
        tmp                 = load(behFile);
        outBeh              = tmp.outBeh;
    else
        
        % process behavioral data
        dt                  = [];
        dt.paths            = paths;
        dt.subject          = subjects{iSub};
        dt.direc            = direc;
        dt.place            = place;
        dt.maxAbsYaw        = param.maxAbsYaw;
        dt.maxNumTrials     = param.maxNumTrials;
        dt.newTimeRes       = param.newTimeRes;
        dt.naviCutoff       = param.naviCutoff;
        dt.naviSmooth       = param.naviSmooth;
        dt.rotSmooth        = param.rotSmooth;
        tmpOut              = LK_EvalBehOverall_20200623(dt);
        
        % use the new timeline for analysis
        outBeh              = tmpOut.new;
        
        % save output
        save(behFile, 'outBeh');
    end
    
    % collect behavioral output for all subjects
    allBeh  = cat(1, allBeh, outBeh);
    
    %% memory performance for each sample
    
    % for each trial, get drop error and normalized memory accuracy
    dropErrorPerTrial       = outBeh.trials(:, end);
    potDropErrorPerTrial    = pdist2(outBeh.trials(:, 9:10), locsInArena);
    memorAccPerTrial        = sum(dropErrorPerTrial < potDropErrorPerTrial, 2) ./ sum(~isnan(potDropErrorPerTrial), 2);
    
    % for each sample, get drop error and normalized memory accuracy
    dropErrorPerSmpl        = nan(size(outBeh.behinfo, 1), 1);
    memorAccPerSmpl         = nan(size(outBeh.behinfo, 1), 1);
    for iTrial = 1:size(dropErrorPerTrial, 1)
        dropErrorPerSmpl(outBeh.behinfo(:, 7) == iTrial)    = dropErrorPerTrial(iTrial, 1);
        memorAccPerSmpl(outBeh.behinfo(:, 7) == iTrial)     = memorAccPerTrial(iTrial, 1);
    end
    
    % sanity check
    if size(unique(memorAccPerTrial(~isnan(memorAccPerTrial))), 1) ~= size(memorAccPerTrial(~isnan(memorAccPerTrial)), 1)
        error('Two or more trials have exactly the same memory accuracy.');
    end
    
    %% mask for analysis timepoints
    
    % boolean that encodes whether specific time bins shall be included in
    % the analysis
    if isfield(param, 'trialPhase') && strcmp(param.trialPhase, 'Active')
        bMask4Analysis  = outBeh.bNavi | outBeh.bRotation;
    end
    
    %% wires to investigate
    
    % microwires
    wires   = dir(fullfile(paths.spike, subjects{iSub}, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    %% loop through wires
    for iWire = 1:size(wires, 1)
        
        % report
        fprintf('\tWire: %s.\n', wires(iWire).name);
        
        %% brain region of this wire
        
        % subject information
        s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion  = s.subjectdata.micro2macro{logIdx, 3};
        
        %% spike times in behavioral time
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
        catch
            fprintf('No wave-clus for this wire.\n');
            exByWC  = cat(1, exByWC, [iSub, iWire]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on inspection)
        c4a     = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', param.clusterName, '.mat'));
        
        % load behavioral times of spikes
        ccb     = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('"cluster_class" not congruent with "cluster_class_behtime".');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
                        
            % skip if this cell was not analyzed in the reference analysis
            if ~any(all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2))
                fprintf('\t\t- This cell was not analyzed in the reference analysis, thus skipping.\n');
                exByRefAna  = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
                continue
            end
            
            % sanity check
            if sum(all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2)) ~= 1
                error('Problem with the reference analysis.');
            end
            
            % data from this cluster
            thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% firing rate for each time bin
            
            % number of spikes and FR for each timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time-bin
            FR          = numSpikes ./ outBeh.durations;
            
            %% empirical evaluation: modulation of FR by memory performance
            
            % reduce FR to relevant time window
            rFR     = FR(bMask4Analysis);
            
            % reduce behavior to relevant time window
            rTime               = outBeh.behinfo(bMask4Analysis, 1); % time
            rMemorAccPerSmpl    = memorAccPerSmpl(bMask4Analysis); % memory accuracy per sample
            
            % if desired, rank memory values
            if strcmp(memor.bRankData, 'yes')
                % rank memory values and normalize maximum to 1
                [~, ~, iM]          = unique(rMemorAccPerSmpl);
                rMemorAccPerSmpl    = iM ./ max(iM);
            end
            
            % correlation
            rhoFR2Mem           = corr(rMemorAccPerSmpl, rFR, 'type', memor.corrType);
            
            % partial correlation to control for potential confounds
            % (time/experience)
            rOtherPred          = rTime; % other predictors
            empRho              = partialcorr(rMemorAccPerSmpl, rFR, rOtherPred, 'type', memor.corrType);
            
            % perform linear regression to get residuals
            mdl                 = fitlm(rOtherPred, rFR);
            
            % estimate empirical values for plotting
            memorLevels         = unique(rMemorAccPerSmpl); % unique memory performances
            empValsRaw          = nan(size(memorLevels, 1), 1); % raw values
            empValsCorr         = nan(size(memorLevels, 1), 1); % corrected values (obtained via residuals from linear regression)
            for iL = 1:size(memorLevels, 1)
                empValsRaw(iL)  = mean(rFR(rMemorAccPerSmpl == memorLevels(iL)));
                empValsCorr(iL) = mean(mdl.Residuals.Raw(rMemorAccPerSmpl == memorLevels(iL))); % effect of other predictors removed
            end
            
            % double-check the number of unique memory levels
            if ~(numel(memorLevels) == numel(unique(outBeh.behinfo(~isnan(outBeh.behinfo(:, 7)) & bMask4Analysis, 7))))
                error('Problem with "memorLevels"');
            end
            
            %% surrogates
            
            % random shifts of surrogate FRs
            randShifts  = datasample(1:numel(rFR) - 1, param.numSurrogates, 'replace', false);
            
            % create surrogates or load previously saved surrogates
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                 = load(resFile);
                rhoFR2MemSurro      = tmp.rhoFR2MemSurro;
                surroRho            = tmp.surroRho;
                surroValsRaw        = tmp.surroValsRaw;
                surroValsCorr       = tmp.surroValsCorr;
            else
                
                % preallocate and loop through surrogates
                rhoFR2MemSurro  = nan(param.numSurrogates, 1); % from correlation
                surroRho        = nan(param.numSurrogates, 1); % from partial correlation
                surroValsRaw    = nan(numel(memorLevels), param.numSurrogates);
                surroValsCorr   = nan(numel(memorLevels), param.numSurrogates);
                parfor iSurro = 1:param.numSurrogates
                    
                    % shift FR in relation to behavior (only include FRs
                    % that were analyzed in the empirical analysis)
                    surroRFR                    = circshift(rFR, randShifts(iSurro));
                    
                    % perform standard correlation
                    rhoFR2MemSurro(iSurro, 1)   = corr(rMemorAccPerSmpl, surroRFR, 'type', memor.corrType);
                    
                    % perform partial correlation
                    surroRho(iSurro, 1)         = partialcorr(rMemorAccPerSmpl, surroRFR, rOtherPred, 'type', memor.corrType);
                    
                    % perform linear regression to get residuals
                    mdl     = fitlm(rOtherPred, surroRFR);
                    
                    % estimate surrogate values
                    tmpRaw  = nan(numel(memorLevels), 1)
                    tmpCorr = nan(numel(memorLevels), 1);
                    for iL = 1:size(memorLevels, 1)
                        tmpRaw(iL, 1)           = mean(surroRFR(rMemorAccPerSmpl == memorLevels(iL)));
                        tmpCorr(iL, 1)          = mean(mdl.Residuals.Raw(rMemorAccPerSmpl == memorLevels(iL)));
                    end
                    surroValsRaw(:, iSurro)     = tmpRaw;
                    surroValsCorr(:, iSurro)    = tmpCorr;
                end
                
                % save surrogates
                save(resFile, ...
                    'rhoFR2MemSurro', 'surroRho', 'surroValsRaw', 'surroValsCorr');
            end
            
            %% interim information
            
            % report rank of empirical value within surrogate values
            fprintf('\n\n==================================================================================\n');
            fprintf('Rank of "empRho" (%.3f) within (i.e., larger than) "surroRho": %.3f.\n', ...
                empRho, sum(empRho > surroRho) / sum(~isnan(surroRho)));
            fprintf('==================================================================================\n\n\n');
            
            %% create figure
            
            % p-value to plot
            if empRho > 0
                pval4plot   = 1 - sum(empRho > surroRho) / sum(~isnan(surroRho)); % positive memory cell?
            else
                pval4plot   = 1 - sum(empRho < surroRho) / sum(~isnan(surroRho)); % negative memory cell?
            end
            
            % figure
            dt                  = [];
            dt.visible          = 'off';
            dt.memorLevels      = memorLevels;
            dt.empValsCorr      = empValsCorr;
            dt.memor.bRankData  = memor.bRankData;
            dt.thisSpike        = thisSpike;
            dt.t.par.sr         = t.par.sr;
            dt.nspk             = sum(numSpikes(bMask4Analysis));
            if pval4plot < 0.01
                dt.figTitle     = {'Memory', ['\itP\rm < 0.01 (', wireRegion, ')']};
            elseif pval4plot < 0.05 && round(pval4plot, 2) == 0.05
                dt.figTitle     = {'Memory', ['\itP\rm < 0.05 (', wireRegion, ')']};
            else
                dt.figTitle     = {'Memory', ['\itP\rm = ', num2str(pval4plot, '%.2f'), ' (', wireRegion, ')']};
            end
            [f, linFit]         = LK_PlotMemoryTuning_20210515(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_FR2Memor'), '-dtiff', '-r450');
                        
            %% collect information for this unit
            
            % basics
            unitRes                             = [];
            unitRes.idx                         = [iSub, iWire, iClus];
            unitRes.wireRegion                  = wireRegion;
                        
            % memory effects
            unitRes.rhoFR2Mem                   = rhoFR2Mem; % empirical correlation
            unitRes.rhoFR2MemSurro              = rhoFR2MemSurro; % surrogate correlations
            unitRes.empRho                      = empRho; % empirical partial correlation
            unitRes.surroRho                    = surroRho; % surrogate partial correlations
            unitRes.memorLevels                 = memorLevels;
            unitRes.empValsRaw                  = empValsRaw;
            unitRes.surroValsRaw                = surroValsRaw;
            unitRes.empValsCorr                 = empValsCorr;
            unitRes.surroValsCorr               = surroValsCorr;
            unitRes.linFit                      = linFit;
            
            % behavioral data
            unitRes.memorAccPerTrial            = memorAccPerTrial;
            
            % detailed data
            unitRes.bMask4Analysis              = bMask4Analysis;
            unitRes.rMemorAccPerSmpl            = rMemorAccPerSmpl;
            unitRes.rFR                         = rFR;
            unitRes.rOtherPred                  = rOtherPred;
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);
            
            %% close all open figures
            close all;
            toc            
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results'));

% unit indices and regions
allUnitIdx  = cell2mat({r.allRes.idx}');
allRegions  = {r.allRes.wireRegion}';

% ELRPD cells from reference analysis
bELPRDCell  = cell2mat({r.ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELPRDCell));
if ~all(all(cell2mat({r.ELRPD.allRes.idx}') == allUnitIdx, 2))
    error('Subject indices from ELRPD-cell analysis different from subject indices of this analysis.\n');
end

% place and direction cells from place/direction-cell analysis
PxDRes      = load('E:\OpenField\PlaceDirAnalysis_030620\20201119_Dir30_Loc10x10\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat');
bDirCell    = cell2mat({PxDRes.allRes.dir_Frank}') > 0.95;
bPlaceCell  = cell2mat({PxDRes.allRes.loc_Frank}') > 0.95;

% report
fprintf('Number of direction cells: %d.\n', sum(bDirCell));
fprintf('Number of place cells: %d.\n', sum(bPlaceCell));

% spatial cells and other, non-spatial cells
bSpatialCell    = bELPRDCell | bDirCell | bPlaceCell;
fprintf('Number of spatial cells: %d.\n', sum(bSpatialCell));
    
%% memory-relevant cells and their overlap with other cell types

fprintf('\nEvaluation of memory-relevant cells.\n');

% collect across cells
allEmpRho       = cell2mat({r.allRes.empRho}'); % from partial correlation
allSurroRho     = cell2mat({r.allRes.surroRho})';

% positive memory cells
posRanks        = sum(allEmpRho > allSurroRho, 2) ./ sum(~isnan(allSurroRho), 2);
bPosMemCell     = posRanks > r.memor.tailThresh;
fprintf('Total number of positive memory cells: %d.\n', sum(bPosMemCell));
fprintf('Number of ELRPD cells that are also positive memory cells: %d.\n', sum(bELPRDCell & bPosMemCell));

% negative memory cells
negRanks        = sum(allEmpRho < allSurroRho, 2) ./ sum(~isnan(allSurroRho), 2);
bNegMemCell     = negRanks > r.memor.tailThresh;
fprintf('Total number of negative memory cells: %d.\n', sum(bNegMemCell));
fprintf('Number of ELRPD cells that are also negative memory cells: %d.\n', sum(bELPRDCell & bNegMemCell));

% memory-relevant cells (either positive or negative)
bMemRelCell     = bPosMemCell | bNegMemCell;
fprintf('Total number of memory-relevant cells: %d.\n', sum(bMemRelCell));

%% figure showing the overlap of memory-relevant cells with spatially-modulated cells

% different experimental groups
groups  = {'ELRPD', 'Memory'; ...
    'Direction', 'Memory'; ...
    'Place-like', 'Memory'; ...
    'Non-spatial', 'Memory'};
nameOfCellType  = {{'Egocentric', 'bearing'}; 'Direction'; 'Place-like'; 'Non-spatial'};
percPerCellType = nan(size(nameOfCellType, 1), 1);

% loop through the different groups and create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [2.9, 2.3, 3.6, 5.6]);
hold on;
for iGroup = 1:size(groups, 1)
        
    % first cell type
    if strcmp(groups{iGroup, 2}, 'Memory')
        bCellA      = bMemRelCell;
        bCellNotA   = ~bMemRelCell;
    end
    
    % second cell type
    if strcmp(groups{iGroup, 1}, 'ELRPD')
        bCellB      = bELPRDCell;
        bCellNotB   = ~bSpatialCell;
    elseif strcmp(groups{iGroup, 1}, 'Direction')
        bCellB      = bDirCell;
        bCellNotB   = ~bSpatialCell;
    elseif strcmp(groups{iGroup, 1}, 'Place-like')
        bCellB      = bPlaceCell;
        bCellNotB   = ~bSpatialCell;
    elseif strcmp(groups{iGroup}, 'Non-spatial')
        bCellB      = ~bSpatialCell;
        bCellNotB   = ~bSpatialCell; % dummy for consistency
    end
    
    % chi-squared test
    n   = [sum(bCellA & bCellB), sum(bCellA & bCellNotB); ...
        sum(bCellNotA & bCellB), sum(bCellNotA & bCellNotB)];
    [X, p]  = myChiSquareTest(n(1, 1), n(1, 2), n(2, 1), n(2, 2));
    fprintf('Chi-squared test to examine interrelation between %s cells (yes vs. other) and memory cells (yes vs. no): X = %.3f, p = %.3f.\n', ...
        groups{iGroup, 1}, X, p);
    disp(n);
    
    % fraction of memory cells in this cell type
    percPerCellType(iGroup, 1)  = sum(bCellA & bCellB) / sum(bCellB);
    
    % bar plot with percentage
    barh(iGroup, 1 * 100, ...
        'FaceColor', [1, 1, 1]); % reference
    barh(iGroup, percPerCellType(iGroup) * 100, ...
        'FaceColor', [0, 0, 0]);
    text(percPerCellType(iGroup) * 100 + 3, iGroup, [num2str(100 * round(percPerCellType(iGroup), 3)), '%'], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, ...
        'verticalalignment', 'middle', 'horizontalalignment', 'left');
end
hold off;
set(gca, ...
    'xlim', [0, 100], ...
    'ylim', [0.4, size(groups, 1) + 0.6], ...
    'ytick', 1:size(groups, 1), 'yticklabel', '', ...
    'ydir', 'reverse', ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
myAxes = get(gca);
for iLabel = 1:size(nameOfCellType, 1)
    text(-8, myAxes.YTick(iLabel), nameOfCellType{iLabel}, 'HorizontalAlignment', 'right', 'fontunits', 'centimeters', 'fontsize', 0.5);
end
xl = xlabel({'% memory cells', 'for each cell type'});
set([gca, xl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'MemRelCellPerc_PerCellType_20210406'), '-dtiff', '-r300');

%% region-wise distribution of memory-relevant cells

fprintf('\nEvaluation of region-wise distribution of memory-relevant cells.\n');

% unique regions and number of cells per unique region
uniqueRegions   = unique(allRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% loop through regions
for iRegion = 1:size(uniqueRegions, 1)
    numCellsPerReg(iRegion, 1)  = sum(strcmp(allRegions, uniqueRegions{iRegion}));
end

% restrict unique regions to those with enough units
uniqueRegions       = uniqueRegions(numCellsPerReg >= 30);

% memory-relevant cells
dt                  = [];
dt.figEx            = [2, 2, 8, 8];
dt.uniqueRegions    = uniqueRegions;
dt.allRegions       = allRegions;
dt.bCell            = bMemRelCell;
dt.ylabel           = 'Memory cells (%)';
[f, percPerReg]     = LK_PlotCellsPerRegion(dt);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'MemRelCell_PerRegion_270120'), '-dtiff', '-r300');
