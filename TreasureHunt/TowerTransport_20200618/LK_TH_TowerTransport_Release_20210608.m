%==========================================================================
% This script examines whether ELRPD cells keep their tuning during
% passive movement.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath('E:\TreasureHunt\Functions\');
addpath(genpath('E:\OpenField\Functions\'));
opengl hardware;

% variables
param               = [];
param.clusterName   = 'Cluster4Analysis_LK_20200416_113912';
param.timeRes       = 0.1; % temporal resolution
param.smoothFac     = 1; % amount of smoothing (1 = no smoothing)
param.numSurrogates = 401; % number of surrogates
param.towerPrePost  = [0, 0]; % additional time before/after
param.corrType      = 'Pearson';
param.myRNG         = 555;
rng(param.myRNG);

% paths
paths       = [];
paths.spike = 'E:\TreasureHunt\SpikeExtraction_NewWC_ManOpt_20200621\';
paths.beh   = 'E:\TreasureHunt\Beh_210420\';
paths.save  = strcat('E:\TreasureHunt\TowerTransport_20200618\20210608\', num2str(1 / param.timeRes), 'Hz_prePost', num2str(unique(param.towerPrePost)), '\');
mkdir(paths.save);

% get subjects
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

%% save settings
save(strcat(paths.save, 'settings'));

%% load previous ELRPD-cell results
ELRPD   = load('E:\TreasureHunt\ELRPDCellAnalysis_220420\20200722_Dir30_Loc6x6_EgoLocDir30\bNaviIfAbove0.1AndbRotationIfAbove0.01_Smooth2\results.mat', ...
    'allRes', 'locDir', 'param');

%% preallocations

% main results for all cells
allRes   	= [];

% bookkeeping
exByWC    	= []; % excluded by wave-clus
exByELRPD  	= []; % excluded by ELRPD cell analysis

%% loop through subjects
for iSub = 1:length(subjects)
    
    % available sessions
    sessions 	= dir(strcat(paths.spike, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        %% behavioral data
        
        % report
        fprintf('\nSubject: %s, session: %s ...\n', subjects{iSub}, sessions(iSess).name);
        
        % original session index
        sessIdx     = split(sessions(iSess).name, '_');
        sessIdx     = str2double(sessIdx{2});
        
        % load behavioral data
        tI          = load(strcat(paths.beh, subjects{iSub}, '\', sessions(iSess).name, '\trialInfo.mat'));
        trialInfo   = tI.trialInfo;
        
        %% wires to investigate
        
        % available microwires
        wires       = dir(fullfile(sessions(iSess).folder, sessions(iSess).name, 'chan*'));
        tmp         = split(transpose({wires.name}), 'n');
        [~, I]      = sort(cellfun(@str2num, tmp(:, 2)));
        wires       = wires(I);
        
        %% loop through wires
        for iWire = 1:size(wires, 1)
            
            % report
            fprintf('\tWire: %s.\n', wires(iWire).name);
            
            %% spike times in behavioral time
            
            % load wave-clus output
            wirePath    = strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep);
            try
                t = load(strcat(wirePath, 'times_datacut.mat'));
            catch
                exByWC  = cat(1, exByWC, [iSub, sessIdx, iWire]); % bookkeeping
                fprintf('No wave-clus for this wire.\n');
                continue;
            end
            
            % load decision whether to use clusters
            c4a = load(strcat(wirePath, param.clusterName, '.mat'));
            
            % load behavioral times of spikes
            ccb = load(strcat(wirePath, 'cluster_class_behtime.mat'));
            
            %% loop through clusters
            for iClus = 1:max(ccb.cluster_class_behtime(:, 1))
                
                % report
                fprintf('\t\tCluster: %d.\n', iClus);
                
                % continue if this unit was not included in the ELRPD-cell
                % analysis
                if ~any(all([iSub, sessIdx, iWire, iClus] == cell2mat({ELRPD.allRes.idx}'), 2))
                    exByELRPD   = cat(1, exByELRPD, [iSub, sessIdx, iWire, iClus]);
                    fprintf('This cluster was not included in the ELRPD-cell analysis.\n');
                    continue;
                end
                
                % data from this cluster
                thisCluster = ccb.cluster_class_behtime(ccb.cluster_class_behtime(:, 1) == iClus, :); % cluster index, microtime (msec), behavioral time (sec)
                
                %% firing rates during tower transport
                
                % time-resolved firing rate during tower transport periods
                onsEnds                 = [trialInfo.TOWER_TRANSPORT_STARTED, trialInfo.TOWER_TRANSPORT_ENDED];
                towerTrialTimeBorders   = cell(size(onsEnds, 1), 1);
                towerTrialFR            = cell(size(onsEnds, 1), 1);
                for iOns = 1:size(onsEnds, 1)
                                        
                    % time-resolved data
                    trialTimeBorders                = (onsEnds(iOns, 1) + param.towerPrePost(1)):param.timeRes:(onsEnds(iOns, 2) + param.towerPrePost(2));
                    towerTrialTimeBorders{iOns, 1}  = trialTimeBorders; % original timeline
                    trialDurs                       = diff(trialTimeBorders, [], 2);
                    trialSpikes                     = histcounts(thisCluster(:, 3), trialTimeBorders);
                    towerTrialFR{iOns, 1}           = trialSpikes ./ trialDurs; % time-resolved firing rate
                end
                
                %% collect information for this unit
                
                % basics
                unitRes                         = [];
                unitRes.idx                     = [iSub, sessIdx, iWire, iClus];
                
                % trial information
                unitRes.trialInfo               = trialInfo;
                
                % tower transport
                unitRes.towerTrialFR            = towerTrialFR;
                unitRes.towerTrialTimeBorders   = towerTrialTimeBorders;
                                                
                % collapse across units
                allRes  = cat(1, allRes, unitRes);                
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));
fprintf('Total number of cells: %d.\n', size(r.allRes, 1));

% indices of all units
allUnitIdx  = cell2mat({r.allRes.idx}');

% ELRPD cells
bELRPDCell  = cell2mat({r.ELRPD.allRes.locdir_clusRank}') > 0.95;
fprintf('Number of ELRPD cells: %d.\n', sum(bELRPDCell));

%% correlation between firing rates and alignments with the preferred ELRPDs during tower transport

% correlation across time, separately for each trial
rho_trFR2prefELRPDAlignment         = cell(size(r.allRes, 1), 1); % empirical correlations
rhoSurro_trFR2prefELRPDAlignment    = cell(size(r.allRes, 1), 1); % surrogate correlations

% loop through cells
for iCell = 1:size(r.allRes, 1)
    
    % process behavior if it is a new subject or a new session
    if iCell == 1 || allUnitIdx(iCell, 1) ~= allUnitIdx(iCell - 1, 1) || allUnitIdx(iCell, 2) ~= allUnitIdx(iCell - 1, 2)
        
        % subject and session index
        subjIdx = allUnitIdx(iCell, 1);
        sessIdx = allUnitIdx(iCell, 2);
        
        % report
        fprintf('Cell: %d, subject: %s, session: %d.\n', iCell, subjects{subjIdx}, sessIdx);
        
        % load behavioral information
        bI     	= load(strcat(paths.beh, subjects{subjIdx}, '\session_', num2str(sessIdx), '\behInfo10Hz.mat'));
        time  	= bI.behInfo.orig.time; % use the data associated with the original time stamps (highest possible resolution)
        xyz    	= bI.behInfo.orig.xyz;
        yaws   	= bI.behInfo.orig.yaws;
    end
    
    % this cell's reference point and preferred ELRPD
    if isempty(r.ELRPD.allRes(iCell).locdir_COMxzClosestLRP)
        continue;
    end
    tmpIdx      = all(r.ELRPD.allRes(iCell).locdir_COMxzClosestLRP == r.ELRPD.locDir.LRPs, 2);
    prefELRPD   = circ_mean(r.ELRPD.locDir.angleCenters', r.ELRPD.allRes(iCell).all_locdir_corrFR(tmpIdx, :)');
    thisCOMxz   = r.ELRPD.allRes(iCell).locdir_COMxz;
    
    % allocentric directions between subject and reference point
    alloDir2COM = atan2(thisCOMxz(1, 2) - xyz(:, 3), thisCOMxz(1, 1) - xyz(:, 1)); % atan2(y, x)
    
    % egocentric directions between subject and reference point
    egoDir2COM 	= angdiff(yaws, alloDir2COM); % angdiff(alpha, beta) subtracts alpha from beta
    
    % alignment of egocentric directions with preferred ELRPD
    prefELRPDAlignment  = cos(angdiff(egoDir2COM, repmat(prefELRPD, size(egoDir2COM, 1), 1)));
        
    % preallocate
    trialRho    	= nan(size(r.allRes(iCell).towerTrialTimeBorders, 1), 1);
    trialRhoSurro 	= nan(size(r.allRes(iCell).towerTrialTimeBorders, 1), r.param.numSurrogates);
    for iOns = 1:size(r.allRes(iCell).towerTrialTimeBorders, 1)
                
        %% trial-wise correlations between firing rates and alignment
        
        % data from this trial
        thisBorders  	= r.allRes(iCell).towerTrialTimeBorders{iOns};
        bThisTrial    	= time >= min(thisBorders) & time < max(thisBorders);
        
        % time-resolved alignment with preferred ELRPD
        trPrefELRPDAlignment    = nan(numel(thisBorders) - 1, 1);
        for iB = 1:numel(thisBorders) - 1
            trPrefELRPDAlignment(iB, 1) = mean(prefELRPDAlignment(time >= thisBorders(iB) & time < thisBorders(iB + 1)));
        end
        
        % correlation between time-resolved alignment and time-resolved
        % firing rates in this trial
        trialRho(iOns, 1)       = corr(trPrefELRPDAlignment, transpose(r.allRes(iCell).towerTrialFR{iOns}), 'type', param.corrType);
        
        % create surrogates by shuffling the firing rates relative to the
        % alignment values
        FR2shuffle              = transpose(r.allRes(iCell).towerTrialFR{iOns});
        [~, randIdx]            = sort(rand(size(FR2shuffle, 1), r.param.numSurrogates));
        shuffledFR              = FR2shuffle(randIdx);
        trialRhoSurro(iOns, :)  = corr(trPrefELRPDAlignment, shuffledFR, 'type', param.corrType);
        
        %% illustration for specific trials
        
        % produce example relationships from trials of ELRPD cells
        if trialRho(iOns, 1) > 0.25 && bELRPDCell(iCell, 1)
            
            %--- relationship between firing rate and alignment with the
            % preferred ELRPD bearing
            f = figure('units', 'centimeters', 'position', [2, 2, 4, 4], 'visible', 'off');
            hold on;
            plot(trPrefELRPDAlignment, r.allRes(iCell).towerTrialFR{iOns}, 'o', ...
                'MarkerEdgeColor', [1, 1, 1], 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerSize', 5);
            P = polyfit(trPrefELRPDAlignment, transpose(r.allRes(iCell).towerTrialFR{iOns}), 1);
            tmpAx = get(gca);
            plot(tmpAx.XLim, P(1) .* tmpAx.XLim + P(2), '-', ...
                'Color', [0, 0, 0], 'LineWidth', 2);
            hold off;
            xl = xlabel('Alignment');
            yl = ylabel('FR (Hz)');
            tl = title(['\itr\rm = ', num2str(trialRho(iOns, 1), '%.2f')]);
            set(gca, ...
                'xlim', tmpAx.XLim, 'xtick', tmpAx.XLim, 'xticklabel', {num2str(tmpAx.XLim(1), '%.1f'), num2str(tmpAx.XLim(2), '%.1f')}, ...
                'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, 'yticklabel', {num2str(tmpAx.YLim(1), '%.0f'), num2str(tmpAx.YLim(2), '%.0f')}, ...
                'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
            set([gca, xl, yl, tl], ...
                'FontUnits', 'centimeters', 'FontSize', 0.5, 'FontWeight', 'normal');
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(r.paths.save, r.subjects{allUnitIdx(iCell, 1)}, '_session_', num2str(allUnitIdx(iCell, 2)), ...
                '_chan', num2str(allUnitIdx(iCell, 3)), '_Clus', num2str(allUnitIdx(iCell, 4)), '_trial', num2str(iOns), ...
                '_trAlignment2towerTrialFR'), '-dtiff', '-r300');
            close(f);
            
            %--- behavioral variables during tower transport
            f = figure('units', 'centimeters', 'position', [2, 2, 8, 8], 'visible', 'off');
            hold on;
            axis square equal;
            % enhance axes
            set(gca, ...
                'xdir', 'reverse', 'ydir', 'reverse', ... % flip x- and y-axis
                'xlim', [min(ELRPD.locDir.xEdges), max(ELRPD.locDir.xEdges)], 'xtick', [min(ELRPD.locDir.xCenters), max(ELRPD.locDir.xCenters)], ...
                'ylim', [min(ELRPD.locDir.zEdges), max(ELRPD.locDir.zEdges)], 'ytick', [min(ELRPD.locDir.zCenters), max(ELRPD.locDir.zCenters)], ...
                'ticklength', [0, 0], 'box', 'on');
            % boundary
            plot(ELRPD.param.arenaCtr(1) + cos(0:0.001:2*pi) .* ELRPD.param.arenaRadius, ELRPD.param.arenaCtr(2) + sin(0:0.001:2*pi) .* ELRPD.param.arenaRadius, '.', ...
                'Color', [0, 0, 0], 'LineWidth', 2);
            
            % yaws at different xz-positions during tower transport
            [qu, qv]    = pol2cart(yaws(bThisTrial), 5 .* ones(size(yaws(bThisTrial)))); % set all radii to 5
            thisXYZ     = xyz(bThisTrial, :);
            lineLength  = 7;
            for iT = 1:size(thisXYZ, 1)
                if mod(iT, 45) == 1 % select timepoints
                    % xz-position
                    plot(thisXYZ(iT, 1), thisXYZ(iT, 3), 'o', ...
                        'Color', [0.7, 0.7, 0.7], 'LineWidth', 1);
                    % allocentric direction towards COM
                    plot([thisXYZ(iT, 1), thisCOMxz(1)], [thisXYZ(iT, 3), thisCOMxz(2)], ':', ...
                        'Color', [0.7, 0.7, 0.7], 'LineWidth', 1);
                    % yaw
                    arrow([thisXYZ(iT, 1), thisXYZ(iT, 3)], [thisXYZ(iT, 1) + qu(iT, 1) .* lineLength, thisXYZ(iT, 3) + qv(iT, 1) .* lineLength], ...
                        10, 'BaseAngle', 50, 'TipAngle', 30, 'width', 1, ...
                        'Color', [0, 0.5, 1]);
                end
            end
            % reference point
            myColors   	= flipud(circshift(hsv, size(hsv, 1) / 2)); % so that red means "ahead" and "right" is associated with green
            potAngles  	= linspace(-pi, pi, size(myColors, 1) + 1);
            prefAngBin  = discretize(prefELRPD, potAngles);
            plot(thisCOMxz(1, 1), thisCOMxz(1, 2), 'o', ...
                'MarkerFaceColor', myColors(prefAngBin, :), 'MarkerEdgeColor', [0, 0, 0], ...
                'MarkerSize', 10, 'LineWidth', 2);
            hold off;
            % improve axes
            xl = xlabel('x (vu)', ...
                'units', 'normalized', 'position', [0.5, -0.025, 0]);
            yl = ylabel('z (vu)', ...
                'units', 'normalized', 'position', [-0.025, 0.5, 0]);
            set([gca, xl, yl], ...
                'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
            ytickangle(gca, 90); % rotate y-tick-labels
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(r.paths.save, r.subjects{allUnitIdx(iCell, 1)}, '_session_', num2str(allUnitIdx(iCell, 2)), ...
                '_chan', num2str(allUnitIdx(iCell, 3)), '_Clus', num2str(allUnitIdx(iCell, 4)), '_trial', num2str(iOns), ...
                '_towerTransportIllustration'), '-dtiff', '-r600');
            
            %--- evolvement of alignment-values over time
            f = figure('units', 'centimeters', 'position', [2, 2, 4.2, 3], 'visible', 'off');
            axes('units', 'centimeters', 'position', [1.41, 0.86, 2.12, 1.86]);
            plot((1:numel(prefELRPDAlignment(bThisTrial))) ./ numel(prefELRPDAlignment(bThisTrial)), prefELRPDAlignment(bThisTrial), '-', ...
                'Color', rgb('darkred'), 'LineWidth', 2);
            set(gca, ...
                'xtick', [0, 1], 'xticklabel', {'Start', 'End'}, ...
                'ylim', [-1, 1], 'ytick', [-1, 0, 1], ...
                'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
            yl = ylabel('Alignment');
            set([gca, yl], ...
                'fontunits', 'centimeters', 'fontsize', 0.5);
            tmpAx = get(gca);
            % save figure
            set(f, 'PaperPositionMode', 'auto');
            print(f, strcat(r.paths.save, r.subjects{allUnitIdx(iCell, 1)}, '_session_', num2str(allUnitIdx(iCell, 2)), ...
                '_chan', num2str(allUnitIdx(iCell, 3)), '_Clus', num2str(allUnitIdx(iCell, 4)), '_trial', num2str(iOns), ...
                '_alignmentOverTime'), '-dtiff', '-r300');
            
            % close all open figures
            close all;
        end
    end
    
    % collect time-resolved correlation across trials
    rho_trFR2prefELRPDAlignment{iCell}      = trialRho;
    rhoSurro_trFR2prefELRPDAlignment{iCell} = trialRhoSurro;
end

%% statistics and figures

% report
fprintf('\nCreating figure for the relationship between preferred ELRPD alignment and firing rates.\n');

% different cell groups
groups  = {'nonELRPDCells', 'ELRPDCells'};
m       = cell(1, size(groups, 2));
ste     = cell(1, size(groups, 2));
dp      = cell(1, size(groups, 2));

% figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 9]);
axes('units', 'centimeters', 'position', [2.25, 2.8 , 5.6, 6]);
hold on;
for iGroup = 1:numel(groups)
    
    % select data
    if strcmp(groups{iGroup}, 'ELRPDCells')
        selData         = cellfun(@nanmean, rho_trFR2prefELRPDAlignment(bELRPDCell)); % average across trials
        selDataSurro    = cell2mat(cellfun(@nanmean, rhoSurro_trFR2prefELRPDAlignment(bELRPDCell), 'uniformoutput', false));
    elseif strcmp(groups{iGroup}, 'nonELRPDCells')
        bIsEmpty        = cellfun(@isempty, rho_trFR2prefELRPDAlignment);
        selData         = cellfun(@nanmean, rho_trFR2prefELRPDAlignment(~bELRPDCell & ~bIsEmpty)); % average across trials
        selDataSurro    = cell2mat(cellfun(@nanmean, rhoSurro_trFR2prefELRPDAlignment(~bELRPDCell & ~bIsEmpty), 'uniformoutput', false));
    end
    
    % statistics
    [~, p, ~, stats]    = ttest(selData);
    fprintf('\nEvaluation for "%s".\n', groups{iGroup});
    fprintf('Are the correlation values between firing rates and alignment with the preferred ELRPD above zero (across cells)? t(%d) = %.3f, p = %.3f.\n', ...
        stats.df, stats.tstat, p);
    empRank             = sum(nanmean(selData) > nanmean(selDataSurro)) / sum(~isnan(nanmean(selDataSurro)));
    fprintf('Are the empirical correlation values higher than the surrogate correlation values? rank = %.3f, p = %.3f.\n', ...
        empRank, 1 - empRank);
    
    % mean and SEM
    m{iGroup}   = nanmean(selData);
    ste{iGroup} = nanstd(selData) / sqrt(sum(~isnan(selData)));
    
    % bar
    bar(iGroup, m{iGroup}, ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.3);
    plot([iGroup, iGroup], [m{iGroup} - ste{iGroup}, m{iGroup} + ste{iGroup}], ...
        'Color', [0, 0, 0], 'LineWidth', 2);
    
    % surrogate distribution
    dp{iGroup}  = nanmean(selDataSurro);
    sp = plot(iGroup + randn(numel(dp{iGroup}), 1) .* 0.05 + 0.2, dp{iGroup}, 'o', ...
        'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [1, 1, 1], ...
        'MarkerSize', 2);
end
set(gca, ...
    'xlim', [0.4, numel(groups) + 0.6], 'xtick', 1:numel(groups), 'xticklabel', '', ...
    'ylim', [-0.02, 0.03], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
myAxes          = get(gca);
myXTickLabel    = {{'Non-', 'egocentric-', 'bearing'}; {'Egocentric', 'bearing'}};
for iLabel = 1:size(myXTickLabel, 1)
    text(myAxes.XTick(iLabel), min(myAxes.YLim) .* 1.06, myXTickLabel{iLabel}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica', 'FontUnits', 'centimeters', 'FontSize', 0.5);
end
xl = xlabel('Cell type', ...
    'units', 'normalized', 'position', [0.5, -0.35]);
yl = ylabel('Pearson''s \itr');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
lg = legend(sp, 'Surrogates', ...
    'location', 'southwest', 'box', 'off');
lg.ItemTokenSize = [3, 18];
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(paths.save, 'ELRPDVSnonELRPDCells_trFR2prefELRPDAlignment_acrossCells'), '-dtiff', '-r450');

% test the two groups against each other
[~, p, ~, stats]   = ttest2(cellfun(@nanmean, rho_trFR2prefELRPDAlignment(~bELRPDCell)), ...
    cellfun(@nanmean, rho_trFR2prefELRPDAlignment(bELRPDCell)));
fprintf('\nAre the correlation values between firing rates and alignment with the preferred ELRPD different between non-ELRPD vs. ELRPD cells? t(%d) = %.3f, p = %.3f.\n', ...
    stats.df, stats.tstat, p);
