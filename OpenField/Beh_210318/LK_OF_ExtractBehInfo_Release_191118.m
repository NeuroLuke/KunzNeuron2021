%==========================================================================
% This script extracts relevant behavioral information from the logfiles of
% the OpenField task.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% path with behavioral data
beh_path = 'E:\OpenField\Beh_210318\';

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
for isub = 1:length(subjects)
    
    % report
    fprintf('\n============================================================\n');
    fprintf('SUBJECT: %s\n', subjects{isub, 1});
    
    % total number of trials
    beh_path_spec   = strcat(beh_path, subjects{isub, 1});
    logfile       	= dir(strcat(beh_path_spec, filesep, '*.log'));
    fid             = fopen(fullfile(logfile.folder, logfile.name));
    C               = textscan(fid, '%s');
    ntrials         = sum(strcmp(C{1}, 'TrialPhase2Counter')); % number of trials
    fclose(fid);
    
    %% read in data
    
    % open logfile
    fid         = fopen(fullfile(logfile.folder, logfile.name));
    textdata  	= [];
    
    % initialize line index
    iline = 1;
    
    % read logfile line-by-line
    tline = fgetl(fid);
    while ischar(tline)
        
        % combine all text-lines into one variable
        textdata{iline, 1} = tline;
        
        % increase line index
        iline = iline + 1;
        
        % read next line
        tline = fgetl(fid);
    end
    
    % close logfile
    fclose(fid);
    
    % report progress
    fprintf('You have read %d text lines...\n', iline);
    
    %% extract trial information ==> "trials.mat"
    
    % preallocate
    trials      = nan(ntrials, 13); % trialidx, object, ITI-onset, cue-onset, retrieval-onset, feedback-onset, reencoding-onset, grab-onset, x-correct, y-correct, x-drop, y-drop, drop error
    trialidx  	= 1;
    
    % loop through lines of textdata and extract information
    for iline = 1:size(textdata, 1)
        
        % report progress
        if mod(iline, 10000) == 0
            fprintf('Percentage looked through: %.3f.\n', 100 * iline / size(textdata, 1));
        end
        
        % information of this specific line
        thislinetext = textdata{iline, :};
        
        % extract information
        if regexp(thislinetext, ' Cue ') % trial index and object-identity
            
            % trial idx
            trials(trialidx, 1)     = trialidx;
            
            % object-identity
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 2)     = str2double(splits{5});
            
        elseif regexp(thislinetext, 'ITI_start') % ITI start
            
            % timepoint of ITI start
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 3)     = str2double(splits{3});

        elseif regexp(thislinetext, 'CUE_start') % cue and retrieval start
            
            % timepoint of cue start and retrieval onset
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 4)     = str2double(splits{3});
            trials(trialidx, 5)     = trials(trialidx, 4) + 2; % cue lasts two seconds until beginning of retrieval
            
        elseif regexp(thislinetext, ' Drop ') % response location
            
            % response location (xy)
            splits                      = strsplit(thislinetext, ' ');
            trials(trialidx, [11, 12])  = [str2double(splits{7}), str2double(splits{9})];
            
        elseif regexp(thislinetext, 'DropBeep') % feedback start
            
            % timepoint of feedback
            splits              = strsplit(thislinetext, ' ');
            trials(trialidx, 6) = str2double(splits{2});
            
        elseif regexp(thislinetext, ' Show ') % correct location
            
            % correct location (xy)
            splits                      = strsplit(thislinetext, ' ');
            trials(trialidx, [9, 10])   = [str2double(splits{7}), str2double(splits{9})];
            
        elseif regexp(thislinetext, 'how accurately placed') % drop error, as estimated by Unreal Engine
            
            % drop error
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 13)    = str2double(splits{6});
            
        elseif regexp(thislinetext, 'STOP_SMILEYFB') % re-encoding onset; cave: sometimes there is no feedback-end before grab
            
            % timepoint of re-encoding start
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 7)     = str2double(splits{4});
            
        elseif regexp(thislinetext, 'GrabBeep') % grab onset
            
            % timepoint of grab
            splits                  = strsplit(thislinetext, ' ');
            trials(trialidx, 8)     = str2double(splits{2});
            
            % if grab is within feedback period
            if isnan(trials(trialidx, 7))
                trials(trialidx, 7) = trials(trialidx, 8); % re-encoding onset identical with grab
            elseif trials(trialidx, 7) < trials(trialidx, 6)
                trials(trialidx, 7) = trials(trialidx, 8); % re-encoding onset identical with grab
            end
            
            % go to next trial
            trialidx                = trialidx + 1;            
        end
    end
    
    %% extract path and heading information ==> "paths.mat" and "headings.mat"
    
    % preallocate
    paths     	= nan(size(textdata, 1), 5); % trialidx, time, x, y, z
    headings   	= nan(size(textdata, 1), 4); % trialidx, time, yaw, orient
    trialidx    = 1;
    
    % loop through lines of textdata
    for iline = 1:size(textdata, 1)
        
        % report progress
        if mod(iline, 10000) == 0
            fprintf('Percentage looked through: %.3f.\n', 100 * iline / size(textdata, 1));
        end
        
        % get information of this specific line
        thislinetext = textdata{iline, :};
        
        % extract information
        if regexp(thislinetext, 'Location X')
            
            % split string into parts
            splits            	= strsplit(thislinetext, ' ');
            paths(iline, :)   	= [trialidx, str2double(splits{3}), str2double(splits{6}), str2double(splits{8}), str2double(splits{10})]; 
        
        elseif regexp(thislinetext, 'YawUnit X')
            
            % split string into parts
            splits           	= strsplit(thislinetext, ' ');
            headings(iline, :)  = [trialidx, str2double(splits{3}), str2double(splits{10}), str2double(splits{12})];
        end        
        
        % go to next trial
        if regexp(thislinetext, 'GrabBeep')
            trialidx        = trialidx + 1;            
        end
    end
    
    % enhance output
    paths     = paths(~isnan(paths(:, 1)), :);
    headings  = headings(~isnan(headings(:, 1)), :);
    if size(paths, 1) ~= size(headings, 1)
        error('Size of "paths" and "headings" does not match.');
    end
    
    % report time period of interest
    fprintf('\tFirst ITI start:     %.3f sec.\n', trials(1, 3));
    fprintf('\tLast ITI end:        %.3f sec.\n', trials(end, 4));
    
    % save it
    save([beh_path_spec, filesep, 'trials'], 'trials');
    save([beh_path_spec, filesep, 'headings'], 'headings');
    save([beh_path_spec, filesep, 'paths'], 'paths');
    
    %% extract comprehensive information for each trial ==> "behinfo.mat"
    
    % preallocate
    behinfo  	= nan(size(textdata, 1), 14); % time, xyz, yaw, orientation, trialidx, trialphase, objectidx, correct xy, drop xy, drop error
    
    % initialize
    trialIdx    = nan;
    trialPhase  = nan;
    objectIdx   = nan;
    correctXY   = nan(1, 2);
    dropXY      = nan(1, 2);
    dropError   = nan;
    
    % loop through text lines
    for iline = 1:size(textdata, 1) - 1
        
        % report progress
        if mod(iline, 10000) == 0
            fprintf('Percentage looked through: %.3f.\n', 100 * iline / size(textdata, 1));
        end
        
        % get information of this and the next line
        thislinetext = textdata{iline, :};
        nextlinetext = textdata{iline + 1, :}; % for yaw and orientation
        
        % update current trial number
        if regexp(thislinetext, ' TrialPhase2Counter ')            
            splits          = strsplit(thislinetext, ' ');
            trialIdx     	= str2double(splits{4});
        end
        
        % update object identity
        if regexp(thislinetext, ' Cue ')            
            splits      	= strsplit(thislinetext, ' ');
            objectIdx       = str2double(splits{5});
        end
        
        % update drop location
        if regexp(thislinetext, ' Drop ')            
            splits          = strsplit(thislinetext, ' ');
            dropXY          = [str2double(splits{7}), str2double(splits{9})];
        end
        
        % update correct location
        if ~isempty(regexp(thislinetext, ' Show ', 'once')) && ~isnan(trialIdx)            
            splits       	= strsplit(thislinetext, ' ');
            correctXY       = [str2double(splits{7}), str2double(splits{9})];
        end
        
        % update drop error
        if regexp(thislinetext, 'how accurately placed')            
            splits       	= strsplit(thislinetext, ' ');
            dropError       = str2double(splits{6});
        end
    
        % update current trial-phase (1 = ITI, 2 = cue, 3 = retrieval, 4 =
        % feedback, 5 = reencoding)
        if regexp(thislinetext, 'ITI_start') % ITI onset            
            trialPhase      = 1;
        elseif regexp(thislinetext, 'CUE_start') % cue onset            
            trialPhase      = 2;            
        elseif regexp(thislinetext, 'bReadyDrop True') % retrieval onset            
            trialPhase      = 3;        
        elseif regexp(thislinetext, 'DropBeep') % feedback onset            
            trialPhase      = 4;            
        elseif regexp(thislinetext, 'STOP_SMILEYFB') % re-encoding onset            
            trialPhase      = 5;            
        elseif regexp(thislinetext, 'GrabBeep') % grab onset            
            trialPhase      = 6;
            
            % reset performance information when trial is over
            correctXY       = nan(1, 2);
            dropXY          = nan(1, 2);
            dropError       = nan;
        end
        
        % get additional information and collect all information whenever
        % subject is at a specific location
        if regexp(thislinetext, 'Location X')
            
            % time and xyz-location information
            splits                  = strsplit(thislinetext, ' ');
            behinfo(iline, 1:4)   	= [str2double(splits{3}), str2double(splits{6}), str2double(splits{8}), str2double(splits{10})];
            
            % heading information (contained in next line)
            splits                  = strsplit(nextlinetext, ' ');
            behinfo(iline, 5:6)     = [str2double(splits{10}), str2double(splits{12})];
            
            % trial-index information
            behinfo(iline, 7)       = trialIdx;
            
            % trial-phase information
            behinfo(iline, 8)       = trialPhase;
            
            % object-index information
            behinfo(iline, 9)       = objectIdx;
            
            % correct-location information
            behinfo(iline, 10:11)   = correctXY;
            
            % drop-location information
            behinfo(iline, 12:13)   = dropXY;
            
            % drop-error information
            behinfo(iline, 14)      = dropError;
        end
    end
    
    % remove lines that do not have a time stamp
    fprintf('Size of "behinfo" before cutting NaNs: \t%d x %d.\n', size(behinfo));
    behinfo = behinfo(~isnan(behinfo(:, 1)), :);
    fprintf('Size of "behinfo" after cutting NaNs: \t%d x %d.\n', size(behinfo));
    
    % save behavioral information
    save([beh_path_spec, filesep, 'behinfo'], 'behinfo');
end

% end report
fprintf('\nExtracted behavioral information for all subjects.\n');
