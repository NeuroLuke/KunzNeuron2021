function [] = lk_trialInfoDiagnostics_110520(cfg, trialInfo)
%
% LK_TRIALINFODIAGNOSTICS_100520 produces several figures to double check
% the extraction of behavioral information from the logfile.
%
% Lukas Kunz, 2021

%% first time and settings

% use the first HOMEBASE_TRANSPORT_STARTED as first time
firstTime   = trialInfo.HOMEBASE_TRANSPORT_STARTED(1);

% general arena settings
xLim        = [310, 430];
zLim        = [300, 420];
arenaCtr    = [mean(xLim), mean(zLim)];
arenaRadius = 50;

%% general: number of chests, number of location recalls, number of object recalls

try    
    % number of chests, location recalls, and object recalls
    numChestsAndRecalls     = [trialInfo.numChestsPerTrial, trialInfo.numLocRecallPerTrial, trialInfo.numObjRecallPerTrial];
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 20]);
    hold on;
    barh(1:size(numChestsAndRecalls, 1), numChestsAndRecalls, ...
        'stacked');
    set(gca, ...
        'xlim', [0, 7], 'xtick', 1:7, ...
        'ydir', 'reverse');
    yl = ylabel('Trial index');
    tl = title({'Number of chests (blue)', ...
        'Number of location recall (red)', ...
        'Number of object recall (yellow)'});
    set([gca, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_numChestsAndRecalls'), '-dtiff', '-r150');
catch
    warning('Could not create figure for number of chests per trial.');
end

%% homebase transport: histogram of durations

try    
    % durations
    durs    = trialInfo.HOMEBASE_TRANSPORT_ENDED - trialInfo.HOMEBASE_TRANSPORT_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Homebase transport: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_homebaseTransport_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for homebase transport durations.');
end

%% navigation: chest locations

try    
    % chest locations
    chestLocs   = trialInfo.chestLoc;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    plot(arenaCtr(1) + cos(0:0.001:2*pi) .* arenaRadius, arenaCtr(2) + sin(0:0.001:2*pi) .* arenaRadius, ...
        'Color', [0, 0, 0]);
    plot(chestLocs(:, 1), chestLocs(:, 3), 'o', ...
        'Color', [1, 1, 1], 'MarkerFaceColor', [0, 0, 1]);
    hold off;
    set(gca, ...
        'xlim', xLim, 'ylim', zLim, ...
        'xdir', 'reverse', 'ydir', 'reverse');
    axis square equal;
    xl = xlabel('x');
    yl = ylabel('z');
    tl = title('Navigation: chest locations');
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_chestLocs'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for chest locations.');
end

%% navigation: treasure labels in German and English

try
    % treasure labels
    allTreasLabel   = cell(size(trialInfo.TREASURE_LABEL_EN, 1), 1);
    for iTreas = 1:size(allTreasLabel, 1)
        allTreasLabel{iTreas, 1}    = [trialInfo.TREASURE_LABEL_GER{iTreas, 1}, ', ', trialInfo.TREASURE_LABEL_EN{iTreas, 1}];
    end
    
    % create figure showing the treasure labels per chest
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 22]);
    axes('units', 'centimeters', 'position', [1, 0.5, 6, 20]);
    text(ones(size(allTreasLabel, 1), 1) .* 0.2, (1:size(allTreasLabel, 1)) ./ size(allTreasLabel, 1), allTreasLabel, ...
        'units', 'data', 'fontunits', 'centimeters', 'fontsize', 0.2);
    set(gca, ...
        'ydir', 'reverse');
    axis off;
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_treasureLabels'), '-dtiff', '-r300');    
catch
    warning('Could not create figure for treasure labels.');
end

%% navigation: durations

try    
    % durations
    durs    = trialInfo.TRIAL_NAVIGATION_ENDED - trialInfo.TRIAL_NAVIGATION_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Trial navigation: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_trialNavigation_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for navigation durations.');
end

%% chest rotation: durations

try
    % durations
    durs    = trialInfo.PLAYER_CHEST_ROTATION_ENDED - trialInfo.PLAYER_CHEST_ROTATION_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Chest rotation: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_chestRotation_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for chest-rotation durations.');
end

%% tower transport: durations

try    
    % durations
    durs    = trialInfo.TOWER_TRANSPORT_ENDED - trialInfo.TOWER_TRANSPORT_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Tower transport: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_towerTransport_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for tower-transport durations.');
end

%% distractor game: durations

try
    % durations and correctness of distractor game response
    durs                = trialInfo.DISTRACTOR_GAME_ENDED - trialInfo.DISTRACTOR_GAME_STARTED;
    correctnessDisGame  = strcmp(trialInfo.correctnessDisGame, 'CorrectBoxAudio');
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Distractor game: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')], ...
        ['Correct: ', num2str(sum(correctnessDisGame)), ' of ', num2str(size(correctnessDisGame, 1))]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_distractorGame_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for distractor-game durations.');
end

%% recall: durations

try    
    % durations
    durs    = trialInfo.RECALL_PHASE_ENDED - trialInfo.RECALL_PHASE_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Recall: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_recallPhase_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for recall durations.');
end

%% recall type

try
    % recall type per trial
    recallType  = trialInfo.recallType;
    numRecalls  = size(recallType, 1);
    
    % create figure showing recall types
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 22]);
    axes('units', 'centimeters', 'position', [1, 0.5, 6, 20]);
    hold on;
    text(ones(numRecalls, 1) .* 0.05, (1:numRecalls) ./ numRecalls, num2str(transpose(1:numRecalls)), ...
        'units', 'data', 'fontunits', 'centimeters', 'fontsize', 0.25);
    text(ones(numRecalls, 1) .* 0.2, (1:numRecalls) ./ numRecalls, recallType, ...
        'units', 'data', 'fontunits', 'centimeters', 'fontsize', 0.25);
    hold off;
    set(gca, ...
        'ydir', 'reverse');
    axis off;
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_recallType'), '-dtiff', '-r300');    
catch
    warning('Could not create figure for recall type.');
end

%% object-cued recall of location: cueing object

try    
    % cueing objects
    cueingObj       = trialInfo.cueingObject;
    numCueingObj    = size(cueingObj, 1);
    
    % create figure showing the cueing objects
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 22]);
    axes('units', 'centimeters', 'position', [1, 0.5, 6, 20]);
    text(ones(numCueingObj, 1) .* 0.2, (1:numCueingObj) ./ numCueingObj, cueingObj, ...
        'units', 'data', 'fontunits', 'centimeters', 'fontsize', 0.25);
    set(gca, ...
        'ydir', 'reverse');
    axis off;
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_locationRecall_cueingObj'), '-dtiff', '-r300');    
catch
    warning('Could not create figure for cueing objects.');
end

%% object-cued recall of location: duration

try    
    % durations
    durs    = trialInfo.timeLocRecall - trialInfo.timeCue4LocRecall;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Location recall: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_locationRecall_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for location-recall durations.');
end

%% object-cued recall of location: drop error

try    
    % time, response and correct locations, drop error
    time            = trialInfo.timeLocRecall;
    responseLocs    = trialInfo.CHOSEN_TEST_POSITION;
    correctLocs     = trialInfo.CORRECT_TEST_POSITION;
    dropError       = sqrt((responseLocs(:, 1) - correctLocs(:, 1)) .^ 2 + (responseLocs(:, 3) - correctLocs(:, 3)) .^ 2);
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    plot(time - firstTime, dropError, '.', ...
        'MarkerSize', 10, 'Color', [0 0 0]);
    myLs = lsline;
    set(myLs, 'Color', [1, 0, 0], 'LineWidth', 2);
    [rho, pval] = corr(time, dropError, 'type', 'Spearman');
    hold off;
    xl = xlabel('Time (sec)');
    yl = ylabel('Drop error (vm)');
    tl = title(['rho = ', num2str(rho, '%.2f'), ', p = ', num2str(pval, '%.3f')]);
    set(gca, ...
        'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_locationRecall_dropError2Time'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for drop error over time.');
end

%% location-cued recall of object: timing

try    
    % trial index, recall onsets and ends
    trialIdx    = trialInfo.trialIdx;
    recallOns   = lk_fillCellArray(trialInfo.OBJECT_RECALL_CHOICE_STARTED, trialInfo.numObjRecallPerTrial);
    recallEnds  = lk_fillCellArray(trialInfo.OBJECT_RECALL_CHOICE_ENDED, trialInfo.numObjRecallPerTrial);
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 16, 16]);
    hold on;
    for iTrial = 1:size(trialIdx, 1)
        thisOns     = recallOns{iTrial};
        thisEnds    = recallEnds{iTrial};
        for iRecall = 1:size(thisOns)
            plot([iTrial, iTrial], [thisOns(iRecall), thisEnds(iRecall)] - firstTime, ...
                'lineWidth', 15);
        end
    end
    hold off;
    xl = xlabel('Trial index');
    yl = ylabel('Time (sec)');
    tl = title(strrep('Object recall: starts and ends', '_', '\_'));
    set([gca, xl, yl, tl], ...
        'FontUnits', 'centimeters', 'FontSize', 0.5, 'FontWeight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_objectRecall_startsAndEnds'), '-dtiff', '-r450');    
catch
    warning('Could not create figure for object-recall timing.');
end

%% location-cued recall of object: recording durations

try    
    % durations
    durs    = trialInfo.RECORDING_ENDED - trialInfo.RECORDING_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Object recall: recording durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_objectRecall_recordingDurations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for recording durations.');
end

%% location-cued recall of object: cueing locations

try    
    % cueing locations and Cortana response
    cueingLocs  = trialInfo.cueingLoc;
    cortanaResp = strcmp(trialInfo.CORTANA_RESPONSE, 'CORRECT');
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
    hold on;
    plot(arenaCtr(1) + cos(0:0.001:2*pi) .* arenaRadius, arenaCtr(2) + sin(0:0.001:2*pi) .* arenaRadius, ...
        'Color', [0, 0, 0]);
    plot(cueingLocs(cortanaResp == 0, 1), cueingLocs(cortanaResp == 0, 3), 'o', ...
        'Color', [1, 1, 1], 'MarkerFaceColor', [1, 0, 0]); % locations leading to incorrect recall
    plot(cueingLocs(cortanaResp == 1, 1), cueingLocs(cortanaResp == 1, 3), 'o', ...
        'Color', [1, 1, 1], 'MarkerFaceColor', [0, 0.5, 0]); % locations leading to correct recall
    hold off;
    set(gca, ...
        'xlim', xLim, 'ylim', zLim, ...
        'xdir', 'reverse', 'ydir', 'reverse');
    axis square equal;
    xl = xlabel('x');
    yl = ylabel('z');
    tl = title('Object recall: cueing locations');
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_objectRecall_cueingLocs'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for cueing locations.');
end

%% temporal retrieval: durations

try    
    % durations
    durs    = trialInfo.TEMPORAL_RETRIEVAL_ENDED - trialInfo.TEMPORAL_RETRIEVAL_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Temporal retrieval: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_temporalRetrieval_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for temporal-retrieval durations.');
end

%% temporal retrieval: options

try    
    % number of options and options
    numOptions  = size(trialInfo.OPTION_A_vs_B, 1);
    options     = cell(numOptions, 1);
    for iTrial = 1:size(options, 1)
        options{iTrial, 1}  = [trialInfo.OPTION_A_vs_B{iTrial, 1}, ', ', trialInfo.OPTION_A_vs_B{iTrial, 2}];
    end
    
    % create figure showing the two options per temporal retrieval
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 22]);
    axes('units', 'centimeters', 'position', [1, 0.5, 6, 20]);
    text(ones(numOptions, 1) .* 0.2, (1:numOptions) ./ numOptions, options, ...
        'units', 'data', 'fontunits', 'centimeters', 'fontsize', 0.25);
    set(gca, ...
        'ydir', 'reverse');
    axis off;
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_temporalRetrieval_options'), '-dtiff', '-r300');    
catch
    warning('Could not create figure for temporal-retrieval options.');
end

%% feedback: durations

try
    % durations
    durs    = trialInfo.FEEDBACK_ENDED - trialInfo.FEEDBACK_STARTED;
    
    % create figure
    f = figure('units', 'centimeters', 'position', [1, 1, 8, 8]);
    histogram(durs);
    set(gca, ...
        'box', 'off');
    xl = xlabel('Duration (sec)');
    yl = ylabel('Count');
    tl = title({'Feedback: durations', ...
        ['Minimum: ', num2str(min(durs), '%.2f')]});
    set([gca, xl, yl, tl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
    % save figure
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(cfg.diagnostics_path, cfg.subject, '_', cfg.session, '_feedback_durations'), '-dtiff', '-r150');    
catch
    warning('Could not create figure for feedback durations.');
end

%% close all open figure
close all;
