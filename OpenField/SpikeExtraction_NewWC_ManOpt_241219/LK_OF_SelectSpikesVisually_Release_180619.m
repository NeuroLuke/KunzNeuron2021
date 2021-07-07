%==========================================================================
% This script determines which units to keep for further analysis based on
% visual inspection.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
opengl hardware;

% paths
spike_path  = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_241219\';

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

% settings
myDate      = '20190618_091853';
saveName  	= ['HardCrit_', myDate];

%% loop through subjects
for isub = 1:length(subjects)
    
    % microwires
    wires   = dir(strcat(spike_path, subjects{isub}, filesep, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    clear tmp;
    
    %% loop through wires
    for iwire = 1:size(wires, 1)
        
        % prellocate
        Cluster4Analysis    = [];
        
        % load spike and cluster data
        try
            t = load(fullfile(wires(iwire).folder, wires(iwire).name, 'times_datacut.mat'));
            c = load(fullfile(wires(iwire).folder, wires(iwire).name, 'cluster_class_behtime.mat'));
            
            % sanity check
            if size(t.cluster_class, 1) ~= size(c.cluster_class_behtime, 1)
                error('"t" and "c" do not have the same number of spikes.');
            end            
        catch
            % save decision
            save(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName), ...
                'Cluster4Analysis'); % save empty structure
            warning('Wave-clus did not extract any spikes.');
            continue;
        end
        
        % load previous decision
        if exist(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'), 'file') > 0
            prevDec = load(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'));
        else
            clear prevDec;
        end
        
        %% loop through clusters
        for iclus = min(t.cluster_class(:, 1)):max(t.cluster_class(:, 1))
            
            % report
            fprintf('\n\n\nSubject: %s. Wire: %s. Cluster: %d.\n', subjects{isub}, wires(iwire).name, iclus);
            
            % general information about this cluster
            tmp         = [];
            tmp.Cluster = iclus; % cluster index
            tmp.Date    = date; % bookkeeping
            
            % select all or specific units for further analysis
            if strcmp(saveName, 'All')
                
                % automatic decision
                if iclus == 0
                    tmp.Decision    = 'no'; % exclude "rest" cluster from further analysis
                else
                    tmp.Decision    = 'yes'; % include all other clusters
                end
            else
                
                % visual decision
                if iclus == 0
                    tmp.Decision    = 'no'; % exclude "rest" cluster from further analysis
                else
                    
                    % decide whether to use this cluster or not
                    if exist('prevDec', 'var')
                        
                        % previous decision
                        tmpPrevDec  = strcmp(prevDec.Cluster4Analysis(cell2mat({prevDec.Cluster4Analysis.Cluster}) == iclus).Decision, 'yes');
                        dec         = input(['Make a decision: shall this cluster be used for analyses (previous decision was ---', num2str(tmpPrevDec), '---)? (1 = yes): ']);
                    else
                        
                        % show wave-clus image
                        f = figure('units', 'normalized', 'position', [0, 0, 1, 1], 'visible', 'off');
                        if iclus <= 3
                            Im  = imread(fullfile(wires(iwire).folder, wires(iwire).name, 'fig2print_datacut.png'));
                        elseif iclus <= 8
                            Im  = imread(fullfile(wires(iwire).folder, wires(iwire).name, 'fig2print_datacuta.png'));
                        else
                            warning('Image missing.')
                        end
                        AxesH = axes('Units', 'pixels', 'position', [1, 1, 1536, 864], 'Visible', 'off');
                        image(Im, 'Parent', AxesH);
                        
                        % make a visual decision
                        dec = input('Make a decision: shall this cluster be used for analyses (previous decision was --- nothing ---)? (1 = yes): ');
                    end
                    if dec == 1
                        tmp.Decision    = 'yes';
                    else
                        tmp.Decision    = 'no';
                    end
                    
                    % close figure
                    close all;                    
                end
            end
            
            % collect across all clusters from this wire
            Cluster4Analysis    = cat(1, Cluster4Analysis, tmp);
        end
        
        % save decision
        save(strcat(wires(iwire).folder, filesep, wires(iwire).name, filesep, 'Cluster4Analysis_', saveName), 'Cluster4Analysis');
    end
end
