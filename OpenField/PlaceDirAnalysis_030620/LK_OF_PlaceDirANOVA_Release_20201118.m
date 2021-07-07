%==========================================================================
% This script detects direction cells and place-like cells via an ANOVA
% procedure.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.myRNG             = 1;
param.clusterName       = 'HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768;
param.newTimeRes        = 0.1;
param.naviCutoff        = 0.001;
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates    	= 101; % number of surrogates
%- for direction
direc                   = [];
direc.angularRes        = 30; % angular resolution
direc.angleEdges        = deg2rad(-180:direc.angularRes:180);
direc.angleCenters      = transpose(movmean(direc.angleEdges, 2, 'endpoints', 'discard'));
%- for location
place                   = [];
place.locRes            = 10; % spatial resolution
place.xEdges            = linspace(-4500, 4500, place.locRes + 1);
place.xCenters          = movmean(place.xEdges, 2, 'endpoints', 'discard');
place.yEdges            = linspace(-4500, 4500, place.locRes + 1);
place.yCenters          = movmean(place.yEdges, 2, 'endpoints', 'discard');
place.cutoffLocBin      = 0.95;
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per bin
ANOVA.model             = 'linear';
ANOVA.sstype            = 2;

% paths
paths       = [];
paths.info  = 'E:\OpenField\SubjectInformation_220318\';
paths.spike = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh   = 'E:\OpenField\Beh_210318\';
paths.save  = strcat('E:\OpenField\PlaceDirAnalysis_030620\20201119_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
    '\bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotation_Smooth', num2str(param.naviSmooth), '\');
mkdir(paths.save);
paths.arena = 'E:\OpenField\Arena\';

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

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByVisInsp     = []; % excluded by visual inspection
exByNumSpikes   = []; % excluded due to low number of spikes

% behavioral information
allBeh          = [];

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
    
    %% behavioral data
    
    % process behavior
    dt              = [];
    dt.paths        = paths;
    dt.subject      = subjects{iSub};
    dt.direc        = direc;
    dt.place        = place;
    dt.maxAbsYaw    = param.maxAbsYaw;
    dt.maxNumTrials = param.maxNumTrials;
    dt.newTimeRes   = param.newTimeRes;
    dt.naviCutoff   = param.naviCutoff;
    dt.naviSmooth   = param.naviSmooth;
    dt.rotSmooth    = param.rotSmooth;
    outBeh          = LK_EvalBehOverall_20200623(dt);
    
    % use new timeline for the analyses
    outBeh          = outBeh.new;
    
    % bookkeeping
    allBeh          = cat(1, allBeh, outBeh);
    
    %% mask for analysis timepoints
    
    % boolean that encodes whether specific time bins shall be included in
    % the analysis
    bMask4Analysis  = outBeh.bNavi | outBeh.bRotation;
    
    %% microwires to investigate
    
    % available microwires
    wires   = dir(fullfile(paths.spike, subjects{iSub}, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    %% loop through wires
    for iWire = 1:size(wires, 1)
        
        % report
        fprintf('\tWire: %s.\n', wires(iWire).name);
        
        %% brain region for this wire
        
        % subject information
        s           = load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx      = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion  = s.subjectdata.micro2macro{logIdx, 3};
        
        %% spike times in behavioral time
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
        catch
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC  = cat(1, exByWC, [iSub, iWire]);
            continue;
        end
        
        % decision whether to use clusters
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', param.clusterName, '.mat'));
        
        % behavioral times of spikes
        ccb = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('Size of original "cluster_class" not congruent with new "cluster_class_behtime".');
        end
        
        %% loop through clusters
        for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
            
            % skip, if necessary
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.Cluster}) == iClus).Decision, 'no')
                fprintf('\t\t- You decided not to analyse this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster     = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime, behavioral time
            thisSpike       = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% number of spikes per time bin
            
            % number of spikes for each behavioral timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time bin
            
            % if the general firing rate is too low, continue
            if (sum(numSpikes(bMask4Analysis)) / sum(outBeh.durations(bMask4Analysis))) < 0.1
                fprintf('The mean firing rate of this cluster is less than 0.1 Hz, thus skipping.\n');
                exByNumSpikes   = cat(1, exByNumSpikes, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            %% ANOVA for effect of direction and place on firing rate
            
            % dependent factor: firing rate
            FR              = numSpikes ./ outBeh.durations;
            
            % compute ANOVAN
            dt              = [];
            dt.FR           = FR(bMask4Analysis);
            dt.group        = [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis)];
            dt.groupNames   = {'D'; 'P'}; % predictors
            dt.ANOVA        = ANOVA;            
            empANOVA        = LK_CompANOVA_070719(dt);
            
            % F values for predictors
            dir_F    = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'D'), strcmp(empANOVA.tbl(1, :), 'F')));
            loc_F    = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'P'), strcmp(empANOVA.tbl(1, :), 'F')));
            
            %% corrected firing rates
            
            % corrected firing rate for direction
            dir_corrFR          = nan(numel(direc.angleCenters), 1);
            tmpFR               = empANOVA.m{strcmp(empANOVA.groupNames, 'D')}(:, 1); % "m" contains estimated means and standard errors
            idx                 = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'D')});
            dir_corrFR(idx)     = tmpFR;
            
            % corrected firing rate for location
            loc_corrFR          = nan(numel(place.yCenters), numel(place.xCenters));
            tmpFR               = empANOVA.m{strcmp(empANOVA.groupNames, 'P')}(:, 1);
            idx                 = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'P')});
            loc_corrFR(idx)     = tmpFR;
            loc_corrFR          = flipud(loc_corrFR); % flip up-down because index 1 shall be in the lower-left corner
            
            %% surrogates
            
            % create surrogates or load previously saved surrogates
            resFile     = strcat(subjSavePath, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved results
                tmp                 = load(resFile);
                dir_Fsurro          = tmp.dir_Fsurro; % surrogate F values for direction
                all_dir_corrFRsurro = tmp.all_dir_corrFRsurro;
                loc_Fsurro          = tmp.loc_Fsurro; % surrogate F value for place
                all_loc_corrFRsurro = tmp.all_loc_corrFRsurro;
            else
                
                % preallocate
                dir_Fsurro          = nan(param.numSurrogates, 1);
                all_dir_corrFRsurro = nan(size(dir_corrFR, 1), param.numSurrogates);
                loc_Fsurro          = nan(param.numSurrogates, 1);
                all_loc_corrFRsurro = nan(size(loc_corrFR, 1), size(loc_corrFR, 2), param.numSurrogates);
                
                % random shifts of firing rates
                randShifts  = transpose(datasample(1:numel(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
                
                % loop through surrogates
                parfor iSurro = 1:param.numSurrogates
                    
                    %% ANOVA
                    
                    % surrogate firing rate
                    FRsurro         = circshift(FR(bMask4Analysis), randShifts(iSurro, 1));
                    
                    % compute ANOVAN
                    dt          	= [];
                    dt.FR         	= FRsurro; % you already masked the FR (see above)
                    dt.group       	= [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis)];
                    dt.groupNames  	= {'D'; 'P'}; % predictors
                    dt.ANOVA       	= ANOVA;
                    surroANOVA    	= LK_CompANOVA_070719(dt);
                    
                    % output F values
                    dir_Fsurro(iSurro, 1)    = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'D'), strcmp(surroANOVA.tbl(1, :), 'F')));
                    loc_Fsurro(iSurro, 1)    = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'P'), strcmp(surroANOVA.tbl(1, :), 'F')));
                    
                    %% corrected firing rates
                    
                    % corrected firing rate for direction
                    dir_corrFRsurro                 = nan(numel(direc.angleCenters), 1);
                    tmpFR                           = surroANOVA.m{strcmp(surroANOVA.groupNames, 'D')}(:, 1);
                    idx                             = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'D')});
                    dir_corrFRsurro(idx)            = tmpFR;
                    all_dir_corrFRsurro(:, iSurro)  = dir_corrFRsurro;
                    
                    % corrected firing rate for location
                    loc_corrFRsurro                     = nan(numel(place.yCenters), numel(place.xCenters));
                    tmpFR                               = surroANOVA.m{strcmp(surroANOVA.groupNames, 'P')}(:, 1);
                    idx                                 = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'P')});
                    loc_corrFRsurro(idx)                = tmpFR;
                    all_loc_corrFRsurro(:, :, iSurro)   = flipud(loc_corrFRsurro); % flip up-down because index 1 shall be in the lower-left corner
                    
                end
                
                % save surrogates
                save(resFile, ...
                    'dir_Fsurro', 'all_dir_corrFRsurro', 'loc_Fsurro', 'all_loc_corrFRsurro');
            end
            
            %% statistical evaluation
            
            % rank of directional F-value
            dir_Frank   = sum(dir_F > dir_Fsurro) / sum(~isnan(dir_Fsurro));
            
            % rank of locational F-value
            loc_Frank   = sum(loc_F > loc_Fsurro) / sum(~isnan(loc_Fsurro));
                        
            % location bins (cf. Ekstrom et al., 2003)
            locBins     = sum(loc_corrFR > all_loc_corrFRsurro, 3) ./ sum(~isnan(all_loc_corrFRsurro), 3) > place.cutoffLocBin;
            
            %% collect information for this unit
            
            % basics
            unitRes              	= [];
            unitRes.idx             = [iSub, iWire, iClus];
            unitRes.totNumSpikes    = sum(numSpikes);
            unitRes.bMask4Analysis  = bMask4Analysis;
            unitRes.wireRegion      = wireRegion;
            
            % direction
            unitRes.dir_F           = dir_F; % F value
            unitRes.dir_Fsurro      = dir_Fsurro; % surrogate F values
            unitRes.dir_Frank       = dir_Frank; % rank of empirical within surrogates
            unitRes.dir_corrFR      = dir_corrFR; % tuning curve
            
            % place
            unitRes.loc_F           = loc_F;
            unitRes.loc_Fsurro      = loc_Fsurro;
            unitRes.loc_Frank       = loc_Frank;
            unitRes.loc_corrFR      = loc_corrFR;
            unitRes.locBins         = locBins;
                        
            % collect across units
            allRes  = cat(1, allRes, unitRes);
                        
            %% figure: direction
            
            % data for figure
            dt                  = [];
            dt.visible          = 'off';
            dt.figEx            = [2, 12, 8.5, 8]; % figure extension
            dt.circleEx         = [2.6, 0.4, 5.5, 5.5];
            dt.waveEx           = [1, 5.5, 2, 2];
            dt.angleCenters     = direc.angleCenters;
            dt.angularRes       = direc.angularRes;
            dt.FR               = dir_corrFR; % tuning curve
            dt.angLabels        = {'E', 'S', 'W', 'N'};
            dt.bReverseY        = true; % whether to flip the y-axis
            dt.thisSpike        = thisSpike;
            dt.t.par.sr       	= t.par.sr;
            dt.nspk          	= sum(numSpikes(bMask4Analysis));
            dt.figTitle         = {'Direction', ...
                ['\itP\rm = ', num2str(1 - dir_Frank, '%.2f'), ' (', wireRegion, ')']};
            if (1 - dir_Frank) < 0.01
                dt.figTitle    	= {'Direction', ...
                    ['\itP\rm < 0.01 (', wireRegion, ')']};
            elseif (1 - dir_Frank) < 0.05 && round(1 - dir_Frank, 2) == 0.05
                dt.figTitle   	= {'Direction', ...
                    ['\itP\rm < 0.05 (', wireRegion, ')']};
            end
            
            % plot figure
            f = LK_PlotFRDir_20200904(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_Dir'), '-dtiff', '-r300');
            
            %% figure: place
            
            % augmentation factor
            augFac          = 20;
            
            % resize x- and y-axes            
            resXCenters     = movmean(linspace(min(place.xEdges), max(place.xEdges), augFac * numel(place.xCenters) + 1), 2, 'endpoints', 'discard');
            resYCenters     = movmean(linspace(min(place.yEdges), max(place.yEdges), augFac * numel(place.yCenters) + 1), 2, 'endpoints', 'discard');
            
            % resize the place map matrix and smooth
            resFR           = myresizem(loc_corrFR, augFac);
            dt              = [];
            dt.FR           = resFR;
            dt.smoothFac    = 21; % smoothing factor
            FRSmooth        = LK_LocSmooth_090719(dt);
            
            % resize location bins
            resLocBins      = myresizem(locBins, augFac);
            
            % data for figure
            dt              = [];
            dt.visible      = 'off';
            dt.figEx        = [2, 5, 8.5, 8]; % figure extension
            dt.circleEx     = [2.6, 0.4, 5.5, 5.5];
            dt.waveEx       = [1, 5.5, 2, 2];
            dt.xCenters     = resXCenters;
            dt.yCenters     = resYCenters;
            dt.FR           = FRSmooth.smFR;
            dt.locBins      = resLocBins;
            dt.path         = outBeh.behinfo(:, 2:3); % navigation path
            dt.thisSpike    = thisSpike;
            dt.t.par.sr     = t.par.sr;
            dt.nspk         = sum(numSpikes(bMask4Analysis));
            dt.figTitle     = {'Place', ...
                ['\itP\rm = ', num2str(1 - loc_Frank, '%.2f'), ' (', wireRegion, ')']};
            if (1 - loc_Frank) < 0.01
                dt.figTitle = {'Place', ...
                    ['\itP\rm < 0.01 (', wireRegion, ')']};
            elseif (1 - loc_Frank) < 0.05 && round(1 - loc_Frank, 2) == 0.05
                dt.figTitle = {'Place', ...
                    ['\itP\rm < 0.05 (', wireRegion, ')']};
            end
            
            % plot figure
            f = LK_PlotFRLoc_20200904(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_Loc'), '-dtiff', '-r300');
            
            %% close all open figures
            close all;
            toc            
        end
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));

% all cell indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\nTotal number of cells: %d.\n', size(allUnitIdx, 1));

%% direction and place cells

% identify direction and place cells
bDirCell    = cell2mat({r.allRes.dir_Frank}') > 0.95;
bPlaceCell  = cell2mat({r.allRes.loc_Frank}') > 0.95;

% report
fprintf('\nNumber of direction cells: %d of %d (P = %.3f).\n', sum(bDirCell), size(bDirCell, 1), ...
    myBinomTest(sum(bDirCell), numel(bDirCell), 0.05));
fprintf('Number of place cells: %d of %d (P = %.3f).\n', sum(bPlaceCell), size(bPlaceCell, 1), ...
    myBinomTest(sum(bPlaceCell), numel(bPlaceCell), 0.05));

% evaluate direction cells
rng(param.myRNG);
dt              = [];
dt.r            = r;
dt.bDirCell     = bDirCell;
dt.bPlaceCell   = bPlaceCell;
LK_EvalDirCells_20210515(dt);

%% region-wise distribution of cells

% all regions, unique regions, and number of cells per unique region
allRegions      = {r.allRes.wireRegion}';
uniqueRegions   = unique(allRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% loop through regions
for iReg = 1:size(uniqueRegions, 1)
    
    % identify cells from this region
    fprintf('\nRegion: %s.\n', uniqueRegions{iReg});
    bThisReg              	= strcmp(allRegions, uniqueRegions{iReg});
    numCellsPerReg(iReg, 1) = sum(bThisReg);
    
    % direction cells
    fprintf('Direction cells: %d (of %d, P = %.3f).\n', ...
        sum(bDirCell & bThisReg), sum(bThisReg), ...
        myBinomTest(sum(bDirCell & bThisReg), sum(bThisReg), 0.05));
    
    % place cells
    fprintf('Place cells: %d (of %d, P = %.3f).\n', ...
        sum(bPlaceCell & bThisReg), sum(bThisReg), ...
        myBinomTest(sum(bPlaceCell & bThisReg), sum(bThisReg), 0.05));
end

% restrict unique regions to those with enough units
uniqueRegions   = uniqueRegions(numCellsPerReg >= 30);

%% bar plot for percentage of cells per region

% direction cells
dt                  = [];
dt.figEx            = [5, 5, 8, 8];
dt.uniqueRegions    = uniqueRegions;
dt.allRegions       = allRegions;
dt.bCell            = bDirCell;
dt.ylabel           = 'Direction cells (%)';
f                   = LK_PlotCellsPerRegion(dt);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'DirCell_PerRegion_090719'), '-dtiff', '-r300');

% place cells
dt                  = [];
dt.figEx            = [5, 5, 8, 8];
dt.uniqueRegions    = uniqueRegions;
dt.allRegions       = allRegions;
dt.bCell            = bPlaceCell;
dt.ylabel           = 'Place-like cells (%)';
f                   = LK_PlotCellsPerRegion(dt);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'PlaceCell_PerRegion_20210409'), '-dtiff', '-r300');
