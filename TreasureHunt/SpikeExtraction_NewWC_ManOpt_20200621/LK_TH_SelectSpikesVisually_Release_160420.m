%==========================================================================
% This script creates a file with the decision whether a given cluster
% shall be used for the analyses based on a visual inspection of its
% wave-clus output.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
opengl hardware;

% paths
spike_path  = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';

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

% settings
myDate      = '20200416_113912';
saveName    = ['LK_', myDate]; % save name

% information about clusters
clusterInfo = [];

%% loop through subjects
for isub = 1:length(subjects)
    
    % report
    fprintf('\nWorking on subject "%s" ...\n', subjects{isub});
    cd(spike_path);
    
    % available sessions
    sess    = dir(strcat(spike_path, subjects{isub}, filesep, 'session*'));
    
    %% loop through sessions
    for isess = 1:size(sess, 1)
        
        % available microwires
        wires   = dir(fullfile(sess(isess).folder, sess(isess).name, 'chan*'));
        tmp     = split(transpose({wires.name}), 'n');
        [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
        wires   = wires(I);
        
        %% loop through wires
        for iwire = 1:size(wires, 1)
            
            % prellocate
            Cluster4Analysis    = [];
            
            % load spike data
            try
                % load wave-clus output
                t   = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_datacut.mat'));
            catch
                % save empty structure if no waveclus output available
                save(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName), 'Cluster4Analysis');
                warning('Wave-clus did not extract any spikes.');
                continue;
            end
            
            % load previous decision
            if exist(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'), 'file') > 0
                prevDec = load(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'));
            else
                clear prevDec; % no previous decision
            end
            
            %% loop through clusters
            for iclus = min(t.cluster_class(:, 1)):max(t.cluster_class(:, 1))
                
                % report
                fprintf('\n\n\nSubject: %s. Session: %s. Wire: %s. Cluster: %d.\n', ...
                    subjects{isub}, sess(isess).name, wires(iwire).name, iclus);
            
                % general information about this cluster
                tmpC4A          = [];
                tmpC4A.Cluster  = iclus; % cluster index
                tmpC4A.Date     = date; % bookkeeping
                
                % if you include all analyzable units
                if strcmp(saveName, 'All')
                    
                    % store decision
                    if iclus == 0
                        tmpC4A.Decision = 'no'; % do not include "rest" cluster
                    else
                        tmpC4A.Decision = 'yes';
                    end                    
                else
                    
                    % perform visual control
                    if iclus == 0
                        tmpC4A.Decision = 'no'; % do not include "rest" cluster
                    else
                        
                        % decide whether to use this cluster or not
                        if exist('prevDec', 'var')
                            
                            % show previous decision
                            tmpPrevDec  = strcmp(prevDec.Cluster4Analysis(cell2mat({prevDec.Cluster4Analysis.Cluster}) == iclus).Decision, 'yes');
                            dec         = input(['Make a decision: shall this cluster be used for analyses (previous decision was ---', ...
                                num2str(tmpPrevDec), '---)? (1 = yes): ']);
                        else
                            
                            % show wave-clus image
                            f = figure('units', 'normalized', 'position', [0, 0, 1, 1]);
                            if iclus <= 3
                                Im  = imread(fullfile(wires(iwire).folder, wires(iwire).name, 'fig2print_datacut.png'));
                            elseif iclus <= 8
                                Im  = imread(fullfile(wires(iwire).folder, wires(iwire).name, 'fig2print_datacuta.png'));
                            else
                                Im  = imread(fullfile(wires(iwire).folder, wires(iwire).name, 'fig2print_datacutb.png'));
                            end
                            AxesH = axes('Units', 'pixels', 'position', [1, 1, 1536, 864], 'Visible', 'off'); % adjust according to screen size
                            image(Im, 'Parent', AxesH);
                            
                            % make a visual decision
                            dec = input('Make a decision: shall this cluster be used for analyses (previous decision was --- nothing ---)? (1 = yes): ');
                        end
                        if dec == 1
                            tmpC4A.Decision = 'yes';
                        else
                            tmpC4A.Decision = 'no';
                        end
                        
                        % close figures
                        close all;
                    end
                end
                
                % collect across all clusters of this wire
                Cluster4Analysis    = cat(1, Cluster4Analysis, tmpC4A);
                
                % collect info about all clusters
                if strcmp(tmpC4A.Decision, 'yes') 
                    dec = 1; % recode decision information
                else
                    dec = 0;
                end
                clusterInfo = cat(1, clusterInfo, [isub, isess, iwire, iclus, dec]);
            end
            
            % save decision
            save(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName), ...
                'Cluster4Analysis');
        end
    end
end

%% save summary output and report

% save
save(strcat(spike_path, saveName, '_clusterInfo'), 'clusterInfo');

% report
fprintf('Number of units offered by wave-clus: %d.\n', size(clusterInfo, 1));
fprintf('Number of units accepted for analysis: %d.\n', sum(clusterInfo(:, 5)));
