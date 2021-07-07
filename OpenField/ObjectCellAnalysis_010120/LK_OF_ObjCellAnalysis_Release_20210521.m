%==========================================================================
% This script identifies object cells.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
opengl hardware;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.myRNG             = 9999;
param.clusterName       = 'HardCrit_20190618_091853';
param.maxNumTrials      = 'all';
param.maxAbsYaw         = 32768; % maximum absolute yaw value in the arena
param.newTimeRes        = 0.1; % temporal resolution
param.naviCutoff        = 0.001; % speed cutoff
param.naviSmooth        = 2;
param.rotSmooth         = 2;
param.numSurrogates     = 101; % number of surrogates
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
place.cutoffLocBin      = 0.95;
%- for object
object                  = [];
object.idx              = transpose(0:7); % possible object labels
object.cutoffObjBin     = 0.95; % cutoff for preferred objects (cf. Ekstrom et al., 2003)
%- for ANOVA procedure
ANOVA                   = [];
ANOVA.minNumObsPerBin   = 5; % minimum number of observations per bin
ANOVA.model             = 'linear'; % type of ANOVA model
ANOVA.sstype            = 2; % type of sum of squares

% paths
paths                   = [];
paths.info              = 'E:\OpenField\SubjectInformation_220318\';
paths.spike             = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.beh               = 'E:\OpenField\Beh_210318\';
paths.arena             = 'E:\OpenField\Arena\';
paths.save              = strcat('E:\OpenField\ObjectCellAnalysis_010120\20200803_Dir', num2str(direc.angularRes), ...
    '_Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), '_Obj\', ...
    'bNaviIfAbove', num2str(param.naviCutoff), 'AndbRotation_Smooth', num2str(param.naviSmooth), '\');
mkdir(paths.save);

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

% load reference results (to include exactly the same units)
ELRPD   = load('E:\OpenField\ELRPDCellAnalysis_20200623\20200617_Dir30_Loc10x10_EgoLocDir30\bNaviIfAbove0.001AndbRotation_Smooth2\results.mat', ...
    'allRes', 'allBeh');

%% save settings
save(strcat(paths.save, 'settings'), '-v7.3');

%% preallocations

% main results for each cell and behavioral data for each subject
allRes       	= [];
allBeh       	= [];

% bookkeeping
exByWC      	= []; % excluded by wave-clus
exByRefAna      = []; % excluded due to low number of spikes

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(param.myRNG);
    
    % subject-specific save path
    subjSavePath    = strcat(paths.save, subjects{iSub}, '\');
    if ~exist(subjSavePath, 'dir')
        mkdir(subjSavePath); % create folder
    end
    
    %% process behavioral data
    
    % process behavior
    behFile     = strcat(subjSavePath, 'outBeh.mat');
    if exist(behFile, 'file')
        
        % load previously processed behavioral data
        fprintf('Loading behavioral data from "%s".\n', behFile);
        tmp                 = load(behFile);
        outBeh              = tmp.outBeh;
    else
        
        % process behavior
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
        outBeh              = LK_EvalBehOverall_20200623(dt);
        
        % use new timeline for the analyses
        outBeh              = outBeh.new;
        
        % save output
        save(behFile, 'outBeh');
    end
    
    % collect behavioral output for all subjects
    allBeh  = cat(1, allBeh, outBeh);
    
    %% mask for analysis timepoints
    
    % boolean that encodes whether specific time bins shall be included in
    % the analysis
    bMask4Analysis  = outBeh.bNavi | outBeh.bRotation;
    fprintf('Information: %.3f%% of the "active time" is included into the analysis.\n', ...
        sum(outBeh.durations(bMask4Analysis)) / sum(outBeh.durations(outBeh.bActive)) * 100);
    
    %% wires to investigate
    
    % available microwires
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
        s         	= load(strcat(paths.info, subjects{iSub}, '\subjectdata_180619.mat'));
        logIdx     	= any(cell2mat(s.subjectdata.micro2macro(:, 1)) == iWire, 2);
        wireRegion 	= s.subjectdata.micro2macro{logIdx, 3};
        
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
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', param.clusterName, '.mat'));
        
        % load behavioral times of spikes
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
                        
            % data for this cluster
            thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime, behavioral time
            thisSpike   = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% exclude cluster if not used in the other analyses
            
            if ~any(all([iSub, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2))
                fprintf('Cluster not included in the ELRPD cell analysis, thus skipping.\n');
                exByRefAna   = cat(1, exByRefAna, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            %% number of spikes and FR per time bin
            
            % number of spikes per timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 3), outBeh.time), 0]); % add 0 spikes for the last time-bin
            
            % FR per timepoint
            FR          = numSpikes ./ outBeh.durations;
            
            %% ANOVAN of main effects: direction, place, and object on firing rate
            
            % compute ANOVAN
            dt              = [];
            dt.FR           = FR(bMask4Analysis); % mask with analysis period
            dt.group        = [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis), outBeh.behinfo(bMask4Analysis, 9)];
            dt.groupNames   = {'D'; 'P'; 'Obj'}; % predictors
            dt.ANOVA        = ANOVA;
            empANOVA        = LK_CompANOVA_070719(dt);
            
            % output F value
            obj_F           = cell2mat(empANOVA.tbl(strcmp(empANOVA.tbl(:, 1), 'Obj'), strcmp(empANOVA.tbl(1, :), 'F')));
            
            %% tuning curve
            
            % corrected firing rate for object
            obj_corrFR              = nan(numel(object.idx), 1);
            tmpFR                   = empANOVA.m{strcmp(empANOVA.groupNames, 'Obj')}(:, 1);
            idx                     = cellfun(@str2num, empANOVA.stats.grpnames{strcmp(empANOVA.groupNames, 'Obj')});
            obj_corrFR(idx + 1)     = tmpFR; % cave: the group names for the factor "object" range between 0 and 7
            obj_corrSEM             = nan(numel(object.idx), 1);
            obj_corrSEM(idx + 1)    = empANOVA.m{strcmp(empANOVA.groupNames, 'Obj')}(:, 2); % standard error of the mean
            
            %% surrogates
            
            % random circular shifts of firing rates
            randShifts  = transpose(datasample(1:length(FR(bMask4Analysis)) - 1, param.numSurrogates, 'replace', false));
            
            % create or load previously saved surrogates
            resFile     = strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_surrogates.mat');
            if exist(resFile, 'file')
                
                % load previously saved surrogates
                tmp                     = load(resFile);
                obj_Fsurro              = tmp.obj_Fsurro; % surrogate F values
                all_obj_corrFRsurro     = tmp.all_obj_corrFRsurro; % surrogate tuning curves
            else
                
                % preallocate
                obj_Fsurro              = nan(param.numSurrogates, 1);
                all_obj_corrFRsurro     = nan(size(obj_corrFR, 1), param.numSurrogates); % e.g., size = 8 x 101                
                parfor iSurro = 1:param.numSurrogates
                    
                    % create surrogate firing rate by circularly shifting
                    % firing rates relative to the behavioral data
                    FRsurro         = circshift(FR(bMask4Analysis), randShifts(iSurro, 1)); % mask with analysis period
                    
                    %% ANOVA: direction, place, object
                    
                    % compute ANOVAN
                    dt              = [];
                    dt.FR           = FRsurro; % you already masked the FR by analysis periods (see above)
                    dt.group        = [outBeh.yawBin(bMask4Analysis), outBeh.xyBin(bMask4Analysis), outBeh.behinfo(bMask4Analysis, 9)];
                    dt.groupNames   = {'D'; 'P'; 'Obj'}; % predictors
                    dt.ANOVA        = ANOVA;
                    surroANOVA      = LK_CompANOVA_070719(dt);
                    
                    % output F value
                    obj_Fsurro(iSurro, 1)   = cell2mat(surroANOVA.tbl(strcmp(surroANOVA.tbl(:, 1), 'Obj'), strcmp(surroANOVA.tbl(1, :), 'F')));
                    
                    %% tuning curve
                    
                    % corrected firing rate for object
                    obj_corrFRsurro                 = nan(numel(object.idx), 1);
                    tmpFR                           = surroANOVA.m{strcmp(surroANOVA.groupNames, 'Obj')}(:, 1);
                    idx                             = cellfun(@str2num, surroANOVA.stats.grpnames{strcmp(surroANOVA.groupNames, 'Obj')});
                    obj_corrFRsurro(idx + 1)        = tmpFR; % cave: the group names for the factor "object" range between 0 and 7
                    all_obj_corrFRsurro(:, iSurro)  = obj_corrFRsurro;
                end
                
                % save surrogates
                save(resFile, ...
                    'obj_Fsurro', 'all_obj_corrFRsurro');
            end
            
            %% statistical evaluation of object tuning
                        
            % rank of object F-value
            obj_Frank   = sum(obj_F > obj_Fsurro) / sum(~isnan(obj_Fsurro));
            
            %% identify preferred objects
            
            % preferred objects
            objBins     = (sum(obj_corrFR > all_obj_corrFRsurro, 2) ./ sum(~isnan(all_obj_corrFRsurro), 2)) > object.cutoffObjBin;
                                    
            %% collect information for this unit
            
            % basics
            unitRes                 = [];
            unitRes.idx             = [iSub, iWire, iClus];
            unitRes.bMask4Analysis  = bMask4Analysis;
            unitRes.wireRegion      = wireRegion;
            
            % main ANOVA: object
            unitRes.obj_F           = obj_F;
            unitRes.obj_Fsurro      = obj_Fsurro;
            unitRes.obj_Frank       = obj_Frank;
            unitRes.obj_corrFR      = obj_corrFR;
            unitRes.objBins         = objBins;
            
            % collect across units
            allRes  = cat(1, allRes, unitRes);
                     
            %% figure: object cells
            
            % data for figure
            dt              = [];
            dt.visible      = 'off';
            dt.figEx        = [2, 5, 8.5, 8];
            dt.boxEx        = [4, 1.35, 4.25, 4.5];
            dt.waveEx       = [1, 5.5, 2, 2];
            dt.object       = object;
            dt.FR           = obj_corrFR; % tuning curve
            dt.SEM          = obj_corrSEM; % SEM of tuning curve
            dt.objBins      = objBins; % preferred objects
            dt.trials       = outBeh.trials;
            dt.thisSpike    = thisSpike; % spike waveforms
            dt.t.par.sr     = t.par.sr; % sampling rate
            dt.nspk         = sum(numSpikes(bMask4Analysis)); % number of spikes
            dt.figTitle     = {'Object tuning', ...
                ['\itP\rm = ', num2str(1 - obj_Frank, '%.2f'), ' (', wireRegion, ')']};
            if (1 - obj_Frank) < 0.01
                dt.figTitle = {'Object tuning', ...
                    ['\itP\rm < 0.01 (', wireRegion, ')']};
            elseif (1 - obj_Frank) < 0.05 && round(1 - obj_Frank, 2) == 0.05
                dt.figTitle = {'Object tuning', ...
                    ['\itP\rm < 0.05 (', wireRegion, ')']};
            end
            
            % plot figure
            f = LK_PlotFRObj_20210515(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(paths.save, subjects{iSub}, '_', wires(iWire).name, '_Clus', num2str(iClus), '_Obj'), '-dtiff', '-r300');
            close(f);
            
            %% close all open figures
            close all;
            toc
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r           = load(strcat(paths.save, 'results.mat'));

% all cells' indices
allUnitIdx  = cell2mat({r.allRes.idx}');
fprintf('\n\n============================================================== Results.\n');
fprintf('Total number of cells: %d.\n', size(allUnitIdx, 1));

%% report number of cells

% object cells
bObjCell          	= cell2mat({r.allRes.obj_Frank}') > 0.95;
numPrefObj        	= sum(cell2mat({r.allRes.objBins})', 2);
bObjCell          	= bObjCell & numPrefObj >= 1; % object cells have at least one preferred object (cf. Qasim et al., 2019)
fprintf('Number of object cells: %d (%.3f%%).\n', sum(bObjCell), 100 * sum(bObjCell) / numel(bObjCell));

% ELRPD cells (from ELRPD cell analysis)
bELRPDCell         	= cell2mat({ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELRPDCell));

% conjunctive ELRPD & object cells
bELRPDandObjCell    = bELRPDCell & bObjCell;
fprintf('Number of ELRPD & object cells: %d.\n', sum(bELRPDandObjCell));

%% results per region

% brain region of each cell
allRegions      = {r.allRes.wireRegion}';

% unique regions
uniqueRegions  	= unique(allRegions);
unitsPerRegion  = nan(numel(uniqueRegions), 1);
for iRegion = 1:numel(uniqueRegions)
    unitsPerRegion(iRegion)    = sum(strcmp(allRegions, uniqueRegions{iRegion}));
end

% results separately for each region
for iRegion = 1:size(uniqueRegions, 1)
    
    % cells in this region
    bThisReg    = strcmp(allRegions, uniqueRegions{iRegion});
    fprintf('\nRegion: %s.\n', uniqueRegions{iRegion});
    fprintf('Number of cells: %d.\n', sum(bThisReg));
    
    % object cells
    fprintf('Number of object cells: %d of %d (P = %.3f).\n', ...
        sum(bObjCell & bThisReg), sum(bThisReg), ...
        myBinomTest(sum(bObjCell & bThisReg), sum(bThisReg), 0.05));
end

% only consider regions with a sufficient number of cells
uniqueRegions   = uniqueRegions(unitsPerRegion >= 30);

%% bar plot for percentage of cells per region

% object cells
dt                  = [];
dt.figEx            = [5, 5, 8, 8];
dt.uniqueRegions    = uniqueRegions;
dt.allRegions       = allRegions;
dt.bCell            = bObjCell;
dt.ylabel           = 'Object cells (%)';
[f, percPerReg]   	= LK_PlotCellsPerRegion(dt);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ObjCell_PerRegion_021019'), '-dtiff', '-r300');

%% evaluate object cells

% characterize object cells
dt           	= [];
dt.r            = r;
dt.bObjCell     = bObjCell;
dt.bELRPDCell   = bELRPDCell;
LK_EvalObjCells_090120(dt);

%% overlap between object cells and ELRPD cells

% chi-squared test
n = [sum(bObjCell & bELRPDCell), sum(bObjCell & ~bELRPDCell); ...
    sum(~bObjCell & bELRPDCell), sum(~bObjCell & ~bELRPDCell)];
[X, p] = myChiSquareTest(n(1, 1), n(1, 2), n(2, 1), n(2, 2));
fprintf('Overlap between object cells and ELRPD cells? X2 = %.3f, p = %.3f.\n', X, p);
disp(n);
