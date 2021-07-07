%==========================================================================
% This script performs spike-detection and -clustering using wave-clus.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\wave_clus-NEW')); % wave-clus v3 (Chaure et al., 2018)

% paths
data_path   = 'E:\OpenField\MicroCut_100518\';
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

%% loop through subjects
parfor isub = 1:size(subjects, 1)
    
    % report
    fprintf('\nWorking on subject %s ...\n', subjects{isub});
    rng(1);
    cd(spike_path);
    
    %% microwires
    
    % wires with trigger and neural data
    trigs   = dir(strcat(data_path, subjects{isub}, filesep, 'ainp*'));
    wires   = dir(strcat(data_path, subjects{isub}, filesep, 'chan*'));
    
    % sort wires
    tmp     = split(transpose({wires.name}), 'n');
    [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    % all files
    files   = [trigs; wires];
    
    % report
    fprintf('\n------------ Performing wave_clus\n');
    
    %% loop through files and extract action potentials
    for ifile = 1:size(files, 1)
        
        % report
        fprintf('\tWorking on file %s ...\n', files(ifile).name);
        
        % specific save path
        spike_path_spec = fullfile(spike_path, subjects{isub}, files(ifile).name, filesep);
        
        % copy previously cut data to specific save path
        mkdir(spike_path_spec);
        copyfile(fullfile(files(ifile).folder, files(ifile).name, filesep, 'datacut.mat'), ...
            spike_path_spec);
        
        %% wave clus
        % needs input file containing "data" and "sr"
        
        % spike sorting and clustering if it is a channel-file
        if ~isempty(regexp(files(ifile).name, 'chan', 'once'))
            
            % get file for wave-clus
            cd(spike_path_spec);
            file2use4wc     = {strcat(spike_path_spec, 'datacut.mat')};
            
            % run wave_clus to extract spikes and do clustering
            try
                % perform spike extraction
                myPar                   = [];
                myPar.detection         = 'neg'; % detection type
                myPar.randomseed        = 1;
                Get_spikes(file2use4wc, 'par', myPar);
                
                % perform spike clustering
                file2use4wc             = {strcat(spike_path_spec, 'datacut_spikes.mat')};
                myPar                   = [];
                myPar.randomseed        = 1;
                myPar.template_sdnum    = 1.5;
                myPar.min_clus          = 60; % minimum cluster size
                myPar.max_clus        	= 10; % maximum cluster number
                myPar.mintemp           = 0.05; % minimum temperature
                Do_clustering(file2use4wc, 'par', myPar);
            catch
                warning('Wave-clus could not be performed.');
            end
            
            %% delete original input-file to save space
            delete(strcat(spike_path_spec, 'datacut.mat'));
        end
    end
end