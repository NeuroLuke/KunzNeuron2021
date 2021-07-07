%==========================================================================
% This script aligns behavioral triggers and microwire triggers to get the
% spike times in behavioral time.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% paths
spike_path  = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
trig_path   = 'E:\TreasureHunt\DataOriginal_010518\';

% subjects
subjects    = {...
    'TH_001'; ...
    'TH_002'; ...
    'TH_003'; ...
    'TH_004'; ...
    'TH_005'; ...
    'TH_006'; ...
    'TH_007'; ...
    'TH_008'; ...
    'TH_009'; ...
    'TH_010'; ...
    'TH_011'; ...
    'TH_012'; ...
    };

%% loop through subjects
for isub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n\n', subjects{isub});
    
    % available sessions
    sess    = dir(strcat(trig_path, subjects{isub}, filesep, 'session*'));
    
    %% loop through sessions
    for isess = 1:size(sess, 1)
        
        % report
        fprintf('\n\nSUBJECT: %s. SESSION: %s.\n\n', subjects{isub}, sess(isess).name);
        
        %% process behavioral trigger information
        
        % read trigger file
        f       = dir(fullfile(sess(isess).folder, sess(isess).name, 'Beh', '*EEGLog.txt'));
        fid     = fopen(fullfile(f.folder, f.name));
        behtrig = textscan(fid, '%f\t%d\t%s');
        
        % behavioral trigger timepoints
        behtrig_timepoints  = behtrig{1}(strcmp(behtrig{3}, 'ON')) ./ 1000; % convert to seconds
        fprintf('Behavioral triggers: t(1) = %.3f sec, t(end) = %.3f sec, range = %.3f sec.\n', behtrig_timepoints(1), behtrig_timepoints(end), range(behtrig_timepoints));
        
        %% process microwire trigger information
        
        % microwire trigger file
        tmpAinp2File    = fullfile(spike_path, subjects{isub}, sess(isess).name, 'ainp2', 'datacut.mat');
        if ~exist(tmpAinp2File, 'file')
            fprintf('\n==== No microwire trigger data for: %s.\n', tmpAinp2File);
            continue;
        end
        
        % microwire trigger data
        tmpAinp2        = load(tmpAinp2File);
        microtrig_data  = transpose(double(tmpAinp2.data));
        microtrig_time  = transpose((0:length(microtrig_data) - 1) ./ tmpAinp2.sr); % start at zero (sec)
                
        % identify microtrigger timepoints
        microtrig_thresh        = 30000; % trigger threshold is 30,000
        microtrig_timepoints    = microtrig_time(diff(microtrig_data < microtrig_thresh) > 0);
        
        % overview figure
        f = figure('units', 'normalized', 'position', [0.05, 0.4, 0.9, 0.4]);
        hold on;
        plot(microtrig_time(1:10:end), microtrig_data(1:10:end), ...
            'Color', rgb('gray'));
        plot(microtrig_timepoints, microtrig_thresh, '*', ...
            'color', rgb('blue'));
        axis tight;
        xlabel('Time (s)');
        ylabel('Voltage');
        title(strrep([subjects{isub}, '. ', sess(isess).name, '. ', 'Microwire trigger channel'], '_', '\_'));
        drawnow;
        
        %% check intertrigger intervals
        
        % inter-trigger-intervals
        microtrig_ITI   = diff(microtrig_timepoints);
        behtrig_ITI     = diff(behtrig_timepoints);
        [rho, pval]     = corr(microtrig_ITI, behtrig_ITI);
        
        if rho < 0.999
            error('Correlation between "microtrig_ITI" and "behtrig_ITI" not sufficiently high.');
        end
        
        % report
        fprintf('Number of microtrigger-ITIs = %d. Number of behtrigger-ITIs = %d.\n', size(microtrig_ITI, 1), size(behtrig_ITI, 1));
        fprintf('Correlation between corresponding ITIs: rho = %.8f, p = %.8f.\n', rho, pval);
        fprintf('The maximum difference between corresponding ITIs = %.6f sec.\n', max(abs(microtrig_ITI - behtrig_ITI)));
        
        %% microwires
        
        % available microwires
        wires   = dir(fullfile(spike_path, subjects{isub}, sess(isess).name, 'chan*'));
        tmpName = split(transpose({wires.name}), 'n');
        [~, I]  = sort(cellfun(@str2num, tmpName(:, 2)));
        wires   = wires(I);
        
        %% loop through wires
        for iwire = 1:size(wires, 1)
            
            % report
            fprintf('\tWire: %s.\n', wires(iwire).name);
            
            % load spike data
            try
                t = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_datacut.mat'));
            catch
                warning('Wave-clus did not extract spikes from this channel.');
                
                % save empty structure as output
                cluster_class_behtime   = [];
                save(fullfile(wires(iwire).folder, wires(iwire).name, filesep, 'cluster_class_behtime'), 'cluster_class_behtime');
                continue;
            end
            
            %% loop through spikes and assign behavioral time
            
            % preallocate: cluster index, microtime (msec), behtime (sec)
            cluster_class_behtime    = nan(size(t.cluster_class, 1), 3);
            for ispike = 1:size(t.cluster_class, 1)
                
                % microwire time of this spike
                thistime    = t.cluster_class(ispike, 2) / 1000; % convert into seconds
                
                % get closest microtriggers
                tmpdiff     = microtrig_timepoints - thistime;
                [~, tmpidx] = max(tmpdiff(tmpdiff < 0));
                
                % continue if spike-time outside of analyzable time
                if isempty(tmpidx) || tmpidx == length(microtrig_timepoints)
                    continue;
                end
                
                % relative timing of spike between the previous and the
                % next trigger
                reltiming   = (thistime - microtrig_timepoints(tmpidx)) / (microtrig_timepoints(tmpidx + 1) - microtrig_timepoints(tmpidx));
                
                % corresponding behavioral timing
                beh_time    = behtrig_timepoints(tmpidx) + reltiming * (behtrig_timepoints(tmpidx + 1) - behtrig_timepoints(tmpidx));
                
                % create output
                cluster_class_behtime(ispike, 1:2)  = t.cluster_class(ispike, :); % cluster-class and microwire time (msec)
                cluster_class_behtime(ispike, 3)    = beh_time; % behavioral time (sec)
            end
            
            % sanity check whether you remove units during processing
            uniqueNames         = unique(cluster_class_behtime(~isnan(cluster_class_behtime(:, 1)), 1));
            controlUniqueNames  = unique(t.cluster_class(~isnan(t.cluster_class(:, 1)), 1));
            if numel(uniqueNames(uniqueNames ~= 0)) ~= numel(controlUniqueNames(controlUniqueNames ~= 0))
                warning('Converting spike times into behavioral times removes --- %d --- units', ...
                    numel(controlUniqueNames(controlUniqueNames ~= 0)) - numel(uniqueNames(uniqueNames ~= 0)));
            end
            
            % save output
            save(fullfile(wires(iwire).folder, wires(iwire).name, 'cluster_class_behtime'), 'cluster_class_behtime');
            
            % report how many spikes are within the analyzable time period
            percSpikesWithinTOI = 100 * sum(~isnan(cluster_class_behtime(:, 1))) / size(cluster_class_behtime, 1);
            fprintf('\t\tBetween first and last trigger: %d spikes (%.2f%% of all spikes).\n', sum(~isnan(cluster_class_behtime(:, 1))), ...
                percSpikesWithinTOI);
        end
        
        % close open figures
        close(f);
    end
end
