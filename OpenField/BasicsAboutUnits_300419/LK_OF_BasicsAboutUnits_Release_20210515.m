%==========================================================================
% This script extracts basic information about the recorded units.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\wave_clus-NEW')); % wave clus version 3 (Chaure et al., 2018)

% settings
param               = [];
param.clusterName   = 'Cluster4Analysis_HardCrit_20190618_091853';
param.peakIdx    	= 20; % peak index of the waveform

% paths
paths         	= [];
paths.raw     	= 'E:\OpenField\MicroCut_100518\';
paths.spike   	= 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
paths.pxd     	= 'E:\OpenField\PlaceDirAnalysis_030620\20201119_Dir30_Loc10x10\bNaviIfAbove0.001AndbRotation_Smooth2\';
paths.save    	= 'E:\OpenField\BasicsAboutUnits_300419\Evaluation_20200902\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

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

% save settings
save(strcat(paths.save, 'settings'));

%% load results from place * direction analysis

% in order to only examine cells that were also included in the analyses
pxdRes          = load(strcat(paths.pxd, 'results.mat'), 'allRes');

%% preallocations

% results for all units
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByNoise       = []; % excluded because noise cluster
exByVisInsp     = []; % excluded based on visual inspection
exByPXD         = []; % excluded from analysis (firing rate too low)

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(999);
    
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
        
        %% get spike times
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
        catch
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC  = cat(1, exByWC, [iSub, iWire]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on visual
        % inspection)
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, param.clusterName, '.mat'));
        
        % load behavioral times of spikes
        ccb = load(fullfile(wires(iWire).folder, wires(iWire).name, 'cluster_class_behtime.mat'));
        
        % sanity check
        if size(t.cluster_class, 1) ~= size(ccb.cluster_class_behtime, 1)
            error('Size of "cluster_class" not congruent with "cluster_class_behtime');
        end
        
        %% create bandpass-filtered raw signal
        
        % load raw data
        raw                 = load(strcat(paths.raw, subjects{iSub}, filesep, wires(iWire).name, '\datacut.mat'));
        
        % filter between 300 and 3,000 Hz using wave-clus
        par                 = [];
        par.sr              = raw.sr;
        par.detect_fmin     = t.par.detect_fmin;
        par.detect_fmax     = t.par.detect_fmax;
        par.detect_order    = t.par.detect_order;
        xf_detect           = spike_detection_filter(raw.data, par);
        
        % calculate STD of the filtered data
        STDnoise            = std(xf_detect);
        
        %% loop through clusters of this wire
        for iClus = 0:max(ccb.cluster_class_behtime(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            
            % continue if this cluster is noise
            if iClus == 0
                exByNoise   = cat(1, exByNoise, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % continue if you decided that this cluster is not sufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.Cluster}) == iClus).Decision, 'no')
                fprintf('\t\t- You decided not to analyse this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % continue if this cluster was not included in the PxD analysis
            if ~any(all([iSub, iWire, iClus] == cell2mat({pxdRes.allRes.idx}'), 2))
                fprintf('\t\t- Was not included in the analysis, thus skipping.\n');
                exByPXD     = cat(1, exByPXD, [iSub, iWire, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster     = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster-number, microtime (msec), behavioral time (sec)
            thisSpike       = t.spikes(ccb.cluster_class_behtime(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% ISI refractoriness
            % = percentage of ISIs < 3ms
            
            ISI                 = diff(thisCluster(:, 2)); % inter-spike-intervals (ms)
            nISI                = size(ISI, 1); % number of ISIs
            percISIlessThan3ms  = 100 * sum(ISI < 3) / nISI; % in percent
            
            %% mean firing rate
            
            meanFR  = size(thisCluster, 1) / (range(thisCluster(:, 2)) / 1000); % (Hz)
            
            %% waveform peak SNR (Faraut et al., 2018)
            % (= ratio between the peak amplitude of the mean waveform and
            % the STD of the noise)            
            
            peakAmpl   	= abs(mean(thisSpike(:, param.peakIdx)));
            peakSNR  	= peakAmpl / STDnoise;            
                        
            %% collect information across units
            
            % this unit's results
            unitRes                     = [];
            unitRes.idx                 = [iSub, iWire, iClus];
            unitRes.percISIlessThan3ms  = percISIlessThan3ms;
            unitRes.meanFR              = meanFR;
            unitRes.peakSNR             = peakSNR;
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);            
        end
    end
end

%% save important output
save(strcat(paths.save, 'results'), '-v7.3');

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));

% unit indices
allUnitIdx  = cell2mat({r.allRes.idx}');

%% units per wire

% number of units per wire
uniqueAllWireIdx    = unique(allUnitIdx(:, 1:2), 'rows', 'stable');
fprintf('Number of wires with at least one unit that was analyzed: %d.\n', ...
    size(uniqueAllWireIdx, 1));
numUnitsPerWire     = nan(size(uniqueAllWireIdx, 1), 1);
for iWire = 1:size(uniqueAllWireIdx, 1)
    numUnitsPerWire(iWire)  = sum(all(uniqueAllWireIdx(iWire, 1:2) == allUnitIdx(:, 1:2), 2));
end

% average number of units per wire
fprintf('Number of units per wire: %.3f +/- %.3f (mean +/- SEM).\n', mean(numUnitsPerWire), ...
    std(numUnitsPerWire) / sqrt(size(numUnitsPerWire, 1)));

% histogram showing the number of units per wire
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(numUnitsPerWire, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'box', 'off', ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Units per wire');
yl = ylabel('Number of wires');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'UnitsPerWire_131019'), '-dtiff', '-r450');

%% ISI refractoriness

% ISI refractoriness from all cells
allPercISILessThan3ms   = cell2mat({r.allRes.percISIlessThan3ms}');

% report
fprintf('Average percentage of ISIs < 3ms across units: %.3f +/- %.3f (mean +/- SEM).\n', ...
    mean(allPercISILessThan3ms), std(allPercISILessThan3ms) / sqrt(size(allPercISILessThan3ms, 1)));
fprintf('Number of units with a percentage of ISIs < 3ms higher than 5%%: %d.\n', ...
    sum(allPercISILessThan3ms >= 5));

% histogram quantifying the ISI refractoriness
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allPercISILessThan3ms, 0:0.5:20, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0, 5], 'xtick', 0:1:20, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('% ISI <3 ms');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ISIRefractoriness_131019'), '-dtiff', '-r450');

%% firing rate

% mean firing rate from all cells
allMeanFR   = cell2mat({r.allRes.meanFR}');

% report
fprintf('Mean firing rate across all units: %.3f +/- %.3f Hz (mean +/- SEM).\n', ...
    mean(allMeanFR), std(allMeanFR) / sqrt(size(allMeanFR, 1)));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allMeanFR, 0:1:25, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xtick', 0:5:25, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Mean FR (Hz)');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'MeanFR_131019'), '-dtiff', '-r450');

%% peak SNR

% peak SNR from all cells
allPeakSNR  = cell2mat({r.allRes.peakSNR}');

% report
fprintf('Average peak SNR: %.3f +/- %.3f.\n', mean(allPeakSNR), ...
    std(allPeakSNR) / sqrt(size(allPeakSNR, 1)));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allPeakSNR, 0:3:33, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Peak SNR');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'PeakSNR_131019'), '-dtiff', '-r450');
