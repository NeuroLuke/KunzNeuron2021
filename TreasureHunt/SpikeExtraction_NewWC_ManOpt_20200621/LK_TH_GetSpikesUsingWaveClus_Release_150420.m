%==========================================================================
% This script performs spike-detection and -clustering using wave-clus.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\wave_clus-NEW')); % wave-clus v3 (Chaure et al., 2018)

% paths
data_path   = 'E:\TreasureHunt\MicroCut_130420\';
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

%% loop through subjects
parfor isub = 1:size(subjects, 1)
    
    % report
    fprintf('\nWorking on subject %s ...\n', subjects{isub});
    rng(1);
    cd(spike_path);
    
    % available sessions
    sess    = dir(strcat(data_path, subjects{isub}, filesep, 'session*'));
    
    %% loop through sessions
    for isess = 1:size(sess, 1)
        
        % report
        fprintf('\nWorking on subject: "%s", session: "%s" ...\n', subjects{isub}, sess(isess).name);
        
        %% microwires
        
        % trigger and neural channels
        trigs   = dir(strcat(sess(isess).folder, filesep, sess(isess).name, filesep, 'ainp*'));
        wires   = dir(strcat(sess(isess).folder, filesep, sess(isess).name, filesep, 'chan*'));
        
        % sort files and combine
        tmp     = split(transpose({wires.name}), 'n');
        [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
        files   = [trigs; wires(I)];
        
        % report
        fprintf('\n------------ Performing wave_clus\n');
        
        %% loop through files and extract action potentials
        for ifile = 1:size(files, 1)
            
            % report
            fprintf('\tWorking on file "%s" ...\n', fullfile(files(ifile).folder, files(ifile).name));
            
            % storage destination
            spike_path_spec = fullfile(spike_path, subjects{isub}, filesep, sess(isess).name, filesep, files(ifile).name, filesep);
            
            % copy original data to storage destination       
            mkdir(spike_path_spec);
            status = copyfile(fullfile(files(ifile).folder, files(ifile).name, filesep, 'datacut.mat'), ...
                spike_path_spec);
            
            %% wave clus
            % needs input file containing "data" and "sr"
            
            if ~isempty(regexp(files(ifile).name, 'chan', 'once')) % if it is a neural channel file
                
                % file for wave-clus
                cd(spike_path_spec);
                file2use4wc = {strcat(spike_path_spec, 'datacut.mat')};
                
                % run wave_clus to extract spikes and do clustering
                try
                    % spike extraction
                    myPar                   = [];
                    myPar.detection         = 'neg'; % detection type
                    myPar.randomseed        = 1;
                    Get_spikes(file2use4wc, 'par', myPar);
                    
                    % spike clustering
                    file2use4wc             = {strcat(spike_path_spec, 'datacut_spikes.mat')};
                    myPar                   = [];
                    myPar.randomseed        = 1;
                    myPar.template_sdnum    = 1.5;
                    myPar.min_clus          = 60; % minimum size of a cluster
                    myPar.max_clus        	= 10; % maximum number of clusters
                    myPar.mintemp           = 0.05; % minimum temperature
                    Do_clustering(file2use4wc, 'par', myPar);
                catch
                    warning('Wave-clus could not be performed.');
                end
                
                %% delete copied original file to save space
                delete(strcat(spike_path_spec, 'datacut.mat'));
            end
        end
    end
end