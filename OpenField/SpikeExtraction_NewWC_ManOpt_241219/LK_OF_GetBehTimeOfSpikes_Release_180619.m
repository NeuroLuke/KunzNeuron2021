%==========================================================================
% This script aligns behavioral triggers and microwire-triggers to get the
% spike times in behavioral time.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% paths
spike_path  = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';
beh_path    = 'E:\OpenField\Beh_210318\';

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

% preallocate for control analyses
allsub_numUnits     = nan(length(subjects), 1); % units per subject
percSpikesWithinTOI = []; % percentage of spikes within time of interest

%% loop through subjects
for isub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n\n', subjects{isub});
    numUnits    = 0; % maximum number of units for further analysis
    
    %% load behavioral and microwire trigger-data
    
    % load behavioral trigger data
    trials                  = load(strcat(beh_path, subjects{isub}, filesep, 'trials.mat'));
    trials                  = trials.trials;
    behtrig_timepoints      = trials(:, 4);
    fprintf('Behavioral triggers: t(1) = %.2f sec, t(end) = %.2f sec, range = %.2f sec.\n', behtrig_timepoints(1), behtrig_timepoints(end), range(behtrig_timepoints));
    
    % load microwire-trigger-data
    tmp                     = load(strcat(spike_path, subjects{isub}, filesep, 'ainp2', filesep, 'datacut.mat'));
    microtrig_data          = transpose(double(tmp.data));
    microtrig_time          = transpose((0:length(microtrig_data) - 1) ./ tmp.sr); % start at zero (sec)
    
    % smooth microtrig-data via convolution to avoid tiny bumps
    microtrig_data          = conv(microtrig_data, ones(11, 1) ./ 11, 'same');
    
    % trigger timepoints
    microtrig_timepoints    = microtrig_time(diff(microtrig_data > 2000) == 1);
    microtrig_timepoints    = microtrig_timepoints(diff([microtrig_timepoints; Inf]) > 4);
    
    % overview figure
    f = figure('units', 'normalized', 'position', [0.05, 0.3, 0.9, 0.5]);
    hold on;
    plot(microtrig_time(1:10:end), microtrig_data(1:10:end));
    plot([microtrig_timepoints, microtrig_timepoints], [min(microtrig_data), max(microtrig_data)], '-', ...
        'color', rgb('orange'));
    axis tight;
    xlabel('Time (sec)');
    ylabel('Voltage (uV)');
    title('Microwire trigger channel');
    
    % double-check whether the number of triggers is identical
    if size(behtrig_timepoints, 1) ~= size(microtrig_timepoints, 1)
        warning('Different number of behavioral and microwire triggers.');
    end
    
    %% get inter-trigger-intervals and find best match
    
    % inter-trigger-intervals
    microtrig_ITI   = diff(microtrig_timepoints);
    behtrig_ITI   	= diff(behtrig_timepoints);
    
    % find best match between ITIs
    xcf = nan(length(microtrig_ITI) - length(behtrig_ITI) + 1, 1);
    for n = 1:length(xcf)
        xcf(n, 1)   = corr(behtrig_ITI, microtrig_ITI(n : n + length(behtrig_ITI) - 1));
    end
    [maxval, maxidx]        = max(xcf);
    microtrig_timepoints    = microtrig_timepoints(maxidx : maxidx + length(behtrig_timepoints) - 1);
    
    % report
    fprintf('When aligned, the ITI-correlation between microtriggers and behtriggers is xcf = %.6f.\n', corr(diff(microtrig_timepoints), diff(behtrig_timepoints)));
    fprintf('There are %d triggers in total.\n', length(microtrig_timepoints));
    
    % add to figure
    plot(microtrig_timepoints, 0.9 .* max(microtrig_data), '*', ...
        'color', rgb('red'));
    hold off;
    drawnow;
    
    % check identified triggers
    input('Are all microwire triggers correctly identified?');
    
    %% spike analysis
    
    % microwires
    wires   = dir(fullfile(spike_path, subjects{isub}, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    %% loop through wires
    for iwire = 1:size(wires, 1)
        
        % report
        fprintf('\tWire: %s.\n', wires(iwire).name);
        
        % load wave-clus output
        try
            t = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_datacut.mat'));
        catch
            warning('Wave-clus did not extract spikes from this channel.\n');
            
            % save empty matrix as output if no data available
            cluster_class_behtime   = [];
            save(fullfile(wires(iwire).folder, wires(iwire).name, 'cluster_class_behtime'), 'cluster_class_behtime');
            continue;
        end
        
        %% loop through spikes and assign behavioral time
        
        % cluster, microtime (msec), behavioral time (sec)
        cluster_class_behtime    = nan(size(t.cluster_class, 1), 3);
        for ispike = 1:size(t.cluster_class, 1)
            
            % microwire time of this spike
            thistime        = t.cluster_class(ispike, 2) / 1000; % (sec)
            
            % closest microtriggers
            tmpdiff         = microtrig_timepoints - thistime;
            [~, tmpidx]     = max(tmpdiff(tmpdiff < 0));
            
            % continue, if spike time not within time of interest (first
            % cue until last cue)
            if isempty(tmpidx) || tmpidx == length(microtrig_timepoints)
                continue;
            end
            
            % relative timing of spike between two consecutive triggers
            reltiming   = (thistime - microtrig_timepoints(tmpidx)) / (microtrig_timepoints(tmpidx + 1) - microtrig_timepoints(tmpidx));
            
            % corresponding behavioral timing
            beh_time    = behtrig_timepoints(tmpidx) + reltiming * (behtrig_timepoints(tmpidx + 1) - behtrig_timepoints(tmpidx));
            
            % collect output
            cluster_class_behtime(ispike, 1:2)  = t.cluster_class(ispike, :); % cluster class and microwire timing
            cluster_class_behtime(ispike, 3)    = beh_time; % behavioral timing
        end
        
        %% relevant output
        
        % save output
        save(fullfile(wires(iwire).folder, wires(iwire).name, 'cluster_class_behtime'), 'cluster_class_behtime');
        
        %% bookkeeping
        
        % count the number of units available for analyses
        uniqueNames = unique(cluster_class_behtime(~isnan(cluster_class_behtime(:, 1)), 1));
        numUnits    = numUnits + numel(uniqueNames(uniqueNames ~= 0)); % ignore "rest" clusters
        
        % percentage of spikes within the time window of interest
        fprintf('\t\tBetween first and last cue: %d spikes (%.2f%% of all spikes).\n', sum(~isnan(cluster_class_behtime(:, 1))), ...
            100 * sum(~isnan(cluster_class_behtime(:, 1))) / size(cluster_class_behtime, 1));        
        percSpikesWithinTOI = cat(1, percSpikesWithinTOI, 100 * sum(~isnan(cluster_class_behtime(:, 1))) / size(cluster_class_behtime, 1));
    end
    
    % collect across sessions
    allsub_numUnits(isub, 1)   = numUnits;
end
