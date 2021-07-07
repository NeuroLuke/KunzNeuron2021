function [] = LK_EvalELRPDCells_20210409(dt)
%
% LK_EvalELRPDCells_20210409 evaluates ELRPD cells.
%
% Input is a structure with fields
%   r           --> results from the ELRPD cell analyis
%   bELRPDCell 	--> logical index to ELRPD cells
%   bDirCell    --> logical index to direction cells
%   direc       --> structure with directional information
%   bPlaceCell  --> logical index to place-like cells
% 
% Lukas Kunz, 2021

% report
fprintf('\n----- Evaluating ELRPD cells.\n');
rng(1);

% font size
myFontSize  = 0.5;

% decide whether to rotate the egocentric directions so that ahead is
% pointing north
rotAngle    = deg2rad(90);

%% number of ELRPD cells for each subject

% report
fprintf('\nEvaluating the number of ELRPD cells per subject and session.\n');

% unit indices
allUnitIdx              = cell2mat({dt.r.allRes.idx}');

% number of ELRPD cells per session
numCellsPerSess         = nan(size(dt.r.subjects, 1), 1); % cells
numELRPDCellsPerSess    = nan(size(dt.r.subjects, 1), 1); % ELRPD cells
for iSess = 1:size(dt.r.subjects, 1)
    numCellsPerSess(iSess, 1)       = sum(allUnitIdx(:, 1) == iSess);
    numELRPDCellsPerSess(iSess, 1)  = sum(allUnitIdx(:, 1) == iSess & dt.bELRPDCell == true);
end

% report
fprintf('Number of sessions with at least one ELRPD cell: %d (out of %d).\n', sum(numELRPDCellsPerSess > 0), numel(numELRPDCellsPerSess));
fprintf('Average number of ELRPD cells per session: %.3f +/- %.3f.\n', ...
    mean(numELRPDCellsPerSess), std(numELRPDCellsPerSess) / sqrt(numel(numELRPDCellsPerSess)));

% number of unique subjects with at least one ELRPD cell
strLength   = cellfun('length', dt.r.subjects);
cutSubjects = cell(size(dt.r.subjects, 1), 1);
for iSub = 1:size(dt.r.subjects, 1)
    if strLength(iSub, 1) > min(strLength)
        cutSubjects{iSub, 1}    = dt.r.subjects{iSub, 1}(1:end-1); % cut final letter
    else
        cutSubjects{iSub, 1}    = dt.r.subjects{iSub, 1};
    end
end
uniqueSubjects          = unique(cutSubjects, 'rows', 'stable');
numCellsPerSubj         = nan(size(uniqueSubjects, 1), 1);
numELRPDCellsPerSubj    = nan(size(uniqueSubjects, 1), 1); % number of ELRPD cells per subject
for iUsub = 1:size(uniqueSubjects, 1)
    numCellsPerSubj(iUsub, 1)       = sum(numCellsPerSess(~cellfun(@isempty, regexp(dt.r.subjects, uniqueSubjects{iUsub}))));
    numELRPDCellsPerSubj(iUsub, 1)  = sum(numELRPDCellsPerSess(~cellfun(@isempty, regexp(dt.r.subjects, uniqueSubjects{iUsub}))));
end
fprintf('Number of unique subjects with at least one ELRPD cell: %d (out of %d).\n', sum(numELRPDCellsPerSubj > 0), numel(numELRPDCellsPerSubj));

%% evaluate distance of COMs (= reference points) to environmental center

% report
fprintf('\nEvaluating the center distance of reference points (i.e., reference-field COMs).\n');

% COMs
allCOMxy                = cell2mat({dt.r.allRes.locdir_COMxy}');
allCOMD2Ctr             = pdist2(allCOMxy, [0, 0]); % center-distances of COMs, all cells
ELRPDCells_allCOMD2Ctr  = allCOMD2Ctr(dt.bELRPDCell, 1); % center-distances of COMs, ELRPD cells only

% step through each center-distance bin and estimate significance
rTotal          = 5000; % radius of entire arena
areaTotal       = pi * rTotal ^ 2; % area of entire arena
rEdges          = 0:250:rTotal;
rCenters        = movmean(rEdges, 2, 'endpoints', 'discard');
numCOMsPerR     = nan(numel(rCenters), 1); % observed counts
numCOMsPerRExp  = nan(numel(rCenters), 1); % expected counts
pCOMsPerR       = nan(numel(rCenters), 1); % p-values from binomial tests
for iR = 1:numel(rCenters)
    
    % observed number of reference points
    innerR                  = rEdges(iR);
    outerR                  = rEdges(iR + 1);
    numCOMsPerR(iR, 1)      = sum(ELRPDCells_allCOMD2Ctr >= innerR & ELRPDCells_allCOMD2Ctr < outerR);
    
    % expected number of reference points
    numCOMsPerRExp(iR, 1)   = numel(ELRPDCells_allCOMD2Ctr) * (pi * outerR ^ 2 - pi * innerR ^ 2) / areaTotal;
    
    % binomial test if empirical value higher than expected value
    if numCOMsPerR(iR, 1) > numCOMsPerRExp(iR, 1)
        pCOMsPerR(iR, 1)    = myBinomTest(numCOMsPerR(iR, 1), sum(dt.bELRPDCell), ...
            numCOMsPerRExp(iR, 1) / sum(dt.bELRPDCell));
    end
end

% Bonferroni correction of p-values
pCOMsPerR   = pCOMsPerR .* numel(rCenters);
idxSig      = find(pCOMsPerR < 0.05);
for iSig = 1:numel(idxSig)
    fprintf('There are significantly more reference points at a distance between %.1f and %.1f vu from the center (P = %.3f).\n', ...
        rEdges(idxSig(iSig)), rEdges(idxSig(iSig) + 1), pCOMsPerR(idxSig(iSig)));
end

% differentiate between center and periphery reference points
borderCtrPeri                   = rTotal / 2;
ELRPDCells_allCOMD2Ctr_Group    = (ELRPDCells_allCOMD2Ctr > borderCtrPeri) + 1; % 1 = close to center; 2 = far from center
fprintf('The border between the two groups (n = %d and n = %d) is at %.1f.\n', ...
    sum(ELRPDCells_allCOMD2Ctr_Group == 1), sum(ELRPDCells_allCOMD2Ctr_Group == 2), borderCtrPeri);

%% figure showing the COM distances to the environmental center

% histogram for distance of COMs to the environmental center
f = figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.55, 1.75, 6, 6]);
hold on;
% histograms of both groups
H1 = histogram(ELRPDCells_allCOMD2Ctr(ELRPDCells_allCOMD2Ctr_Group == 1), rEdges, ...
    'facecolor', rgb('darkgreen'), 'edgecolor', 'none', 'facealpha', 1); % center group
H2 = histogram(ELRPDCells_allCOMD2Ctr(ELRPDCells_allCOMD2Ctr_Group == 2), rEdges, ...
    'facecolor', rgb('limegreen'), 'edgecolor', 'none', 'facealpha', 1); % periphery group
H3 = stairs(rEdges, [numCOMsPerRExp; numCOMsPerRExp(end)], ...
    'Color', [0.5, 0.5, 0.5]); % expected distribution
axis square;
set(gca, ...
    'xlim', [min(rEdges) - 100, max(rEdges) + 100], 'xtick', [0, 2500, 5000], ...
    'ylim', [0, max([H1.BinCounts, H2.BinCounts, H3.YData]) + 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
tmpAx = get(gca);
% indicate significance
LK_SigLine(rCenters, [max(tmpAx.YLim); max(tmpAx.YLim) - range(tmpAx.YLim) * 0.025], pCOMsPerR < 0.05); 
hold off;
xl = xlabel('Distance to center (vu)');
yl = ylabel('Count');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5)
uistack(gca, 'top');
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'ELRPDCells_COMxyDistance2Center_120919'), '-dtiff', '-r300');

%% evaluate spatial distribution of reference points

% report
fprintf('\nEvaluating the spatial distribution of reference points (i.e., reference-field COMs).\n');

% reference points of ELRPD cells
ELRPDCells_allCOMxy     = allCOMxy(dt.bELRPDCell, :);

% plot centers of mass within the environment
f = figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
hold on;
plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 5);
plot(cos(0:0.001:(2*pi)) .* borderCtrPeri, sin(0:0.001:(2*pi)) .* borderCtrPeri, ':', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% reference points
scatter(ELRPDCells_allCOMxy(ELRPDCells_allCOMD2Ctr_Group == 1, 1), ...
    ELRPDCells_allCOMxy(ELRPDCells_allCOMD2Ctr_Group == 1, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % center reference points
scatter(ELRPDCells_allCOMxy(ELRPDCells_allCOMD2Ctr_Group == 2, 1), ...
    ELRPDCells_allCOMxy(ELRPDCells_allCOMD2Ctr_Group == 2, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % periphery reference points
hold off;
set(gca, ...
    'ydir', 'reverse', ... % negative y-values at top
    'xlim', [min(dt.r.locDir.xEdges), max(dt.r.locDir.xEdges)], 'xtick', [min(dt.r.locDir.xCenters), max(dt.r.locDir.xCenters)], ...
    'ylim', [min(dt.r.locDir.yEdges), max(dt.r.locDir.yEdges)], 'ytick', [min(dt.r.locDir.yCenters), max(dt.r.locDir.yCenters)], ...
    'tickdir', 'out', 'ticklength', [0 0], 'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
ytickangle(90);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'ELRPDCells_allCOMxy_120919'), '-dtiff', '-r300');

%% session-specific depiction of reference points

% loop through sessions
for iSess = min(allUnitIdx(:, 1)):max(allUnitIdx(:, 1))
    
    % mask for reference cells from this session
    bMask   = dt.bELRPDCell & allUnitIdx(:, 1) == iSess;
        
    % plot centers of mass within the environment
    f = figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
    axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
    hold on;
    % boundary and separation between center and periphery
    plot(cos(0:0.001:(2*pi)) .* 5000, sin(0:0.001:(2*pi)) .* 5000, '-', ...
        'Color', [0, 0, 0], 'LineWidth', 5);
    plot(cos(0:0.001:(2*pi)) .* borderCtrPeri, sin(0:0.001:(2*pi)) .* borderCtrPeri, ':', ...
        'Color', [0, 0, 0], 'LineWidth', 2);
    % reference points
    scatter(allCOMxy(bMask & allCOMD2Ctr < borderCtrPeri, 1), ...
        allCOMxy(bMask & allCOMD2Ctr < borderCtrPeri, 2), 50, 'o', ...
        'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % center COMs
    scatter(allCOMxy(bMask & allCOMD2Ctr > borderCtrPeri, 1), ...
        allCOMxy(bMask & allCOMD2Ctr > borderCtrPeri, 2), 50, 'o', ...
        'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1); % periphery COMs
    hold off;
    set(gca, ...
        'ydir', 'reverse', ... % negative y-values at top
        'xlim', [min(dt.r.locDir.xEdges), max(dt.r.locDir.xEdges)], 'xtick', [min(dt.r.locDir.xCenters), max(dt.r.locDir.xCenters)], ...
        'ylim', [min(dt.r.locDir.yEdges), max(dt.r.locDir.yEdges)], 'ytick', [min(dt.r.locDir.yCenters), max(dt.r.locDir.yCenters)], ...
        'tickdir', 'out', 'ticklength', [0 0], 'box', 'on');
    xl = xlabel('x (vu)', ...
        'units', 'normalized', 'position', [0.5, -0.025, 0]);
    yl = ylabel('y (vu)', ...
        'units', 'normalized', 'position', [-0.025, 0.5, 0]);
    ytickangle(90);
    set([gca, xl, yl], ...
        'fontunits', 'centimeters', 'fontsize', 0.5);
    set(f, 'PaperPositionMode', 'auto');
    print(f, strcat(dt.r.paths.save, dt.r.subjects{iSess}, '_ELRPDCells_allCOMxy_20200812'), '-dtiff', '-r300');
end

%% are reference points from the same session significantly close to each other?

% reset rng
rng(dt.r.param.myRNG);

% number of surrugates
thisNumSurrogates   = 1001;

% reference points of ELRPD cells and their associated session indices
ELRPDCells_allCOMxy = allCOMxy(dt.bELRPDCell, :);
ELRPDCells_sessIdx  = allUnitIdx(dt.bELRPDCell, 1);

% in each session, estimate the average distance between the reference
% points
Demp    = nan(max(ELRPDCells_sessIdx), 1); % empirical distance
for iSess = 1:max(ELRPDCells_sessIdx)

    % distances between all reference points in this session
    tmpD    = pdist2(ELRPDCells_allCOMxy(ELRPDCells_sessIdx == iSess, :), ELRPDCells_allCOMxy(ELRPDCells_sessIdx == iSess, :));
    bMask  	= tril(tmpD) ~= 0;
    
    % average distance in this session
    Demp(iSess, 1)  = mean(tmpD(bMask));
end

% create surrogates by randomly assigning session indices to reference points
Dsurro  = nan(max(ELRPDCells_sessIdx), thisNumSurrogates); % surrogate distances
for iSurro = 1:thisNumSurrogates
    
    % shuffle the session indices relative to the reference points
    shuffSessIdx    = datasample(ELRPDCells_sessIdx, numel(ELRPDCells_sessIdx), 'replace', false);
    
    % in each surrogate session, estimate the average distance between the
    % reference points
    for iSess = 1:max(shuffSessIdx)
        
        % distances between all reference points in this surrogate session
        tmpD  	= pdist2(ELRPDCells_allCOMxy(shuffSessIdx == iSess, :), ELRPDCells_allCOMxy(shuffSessIdx == iSess, :));
        bMask  	= tril(tmpD) ~= 0;
        
        % average distance in this surrogate session
        Dsurro(iSess, iSurro)   = mean(tmpD(bMask));
    end
end

% determine rank of the empirical distance within surrogate distances
tmpRank = sum(nanmean(Demp, 1) < nanmean(Dsurro, 1), 2) ./ sum(~isnan(nanmean(Dsurro, 1)), 2);
fprintf('Is the empirical distance between the reference points from the same session smaller than expected by chance? rank = %.3f, p = %.3f.\n', tmpRank, 1 - tmpRank);

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
hold on;
h1 = histogram(nanmean(Dsurro, 1), ...
    'EdgeColor', 'none', 'FaceColor', [0.7, 0.7, 0.7], ...
    'normalization', 'probability');
plot([nanmean(Demp, 1), nanmean(Demp, 1)], [0, max(h1.Values)], '-', ...
    'Color', [1, 0, 0], 'LineWidth', 2);
xl = xlabel({'Distance between', 'reference points (vu)'});
yl = ylabel('Probability');
set(gca, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'DistanceBetweenRefPointsFromSameSession_20200825'), '-dtiff', '-r300');

%% distribution of preferred ELRPDs (evaluation via circular mean)

% for each cell, get preferred ELRPD
prefELRPD   = nan(size(dt.r.allRes, 1), 1);
for iCell = 1:size(prefELRPD, 1)
    % if the cell has a COM, extract preferred ELRPD
    if ~isempty(dt.r.allRes(iCell).locdir_COMxyClosestLRP)
        thisCOMLRP              = all(dt.r.locDir.LRPs == dt.r.allRes(iCell).locdir_COMxyClosestLRP, 2);
        prefELRPD(iCell, 1)     = circ_mean(dt.r.locDir.angleCenters, ...
            transpose(dt.r.allRes(iCell).locdir_corrFR(thisCOMLRP == true, :)));
    end
end

% examine different groups of ELRPD cells and different angular spaces
groups      = {'ELRPDCellsCtrCOM'; 'ELRPDCellsPeriCOM'};
degSpace    = [360; 180]; % different angular spaces

% loop through groups
for iGroup = 1:size(groups, 1)
    
    % group of cells to analyze
    if strcmp(groups{iGroup}, 'ELRPDCellsCtrCOM')
        
        % ELRPD cells with center reference points
        bCellMask	= dt.bELRPDCell & allCOMD2Ctr < borderCtrPeri;
        
    elseif strcmp(groups{iGroup}, 'ELRPDCellsPeriCOM')
        
        % ELRPD cells with periphery reference points
        bCellMask   = dt.bELRPDCell & allCOMD2Ctr > borderCtrPeri;
    end
    
    % loop through different angular spaces
    for iD = 1:size(degSpace, 1)
                
        % Rayleigh test in specific angular space
        [pval, z] = circ_rtest(prefELRPD(bCellMask) .* (360 / degSpace(iD)));
        fprintf('Rayleigh test for group "%s" in %d°-space: z = %.3f, P = %.3f.\n', ...
            groups{iGroup}, degSpace(iD), z, pval);
        
        % extract number of observations per angular bin
        tmpF    = figure('visible', 'off');
        myP     = polarhistogram(prefELRPD(bCellMask) .* (360 / degSpace(iD)), dt.r.locDir.angleEdges);
        NcircM  = myP.Values;
        close(tmpF);
        
        %--- create figure
        f = figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
        axis off;
        
        % angle-circle
        ax1 = axes('units', 'centimeters', 'position', [1.5, 1.5, 5.5, 5.5]);
        hold on;
        % plot angle lines
        for iAng = 0:30:165
            plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
                'Color', rgb('gray'), 'LineWidth', 1);
        end
        % plot circle
        tmpx = (-pi:0.001:pi) + rotAngle; % rotate so that 0° is pointing north
        tmpc = linspace(0, 1, length(tmpx));
        patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
            [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
            'EdgeColor', 'interp');
        % plot angle text
        t1 = text(1.125, 0, 'L'); % L here, because x-axis will be reversed
        t2 = text(0, 1.15, 'A');
        t3 = text(-1.125, 0, 'R'); % R here, because x-axis will be reversed
        t4 = text(0, -1.15, 'B');
        set([t1, t2, t3, t4], ...
            'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
            'HorizontalAlignment', 'center');
        % indicate maximum
        text(0.175, 0.15, ['max ', num2str(max(NcircM), '%.0f')], ...
            'units', 'normalized', ...
            'horizontalalignment', 'right', 'verticalalignment', 'top', ...
            'fontunits', 'centimeters', 'fontsize', myFontSize)
        hold off;
        % cave: reverse x-axis so that -90° (i.e., left) points to the left
        set(gca, ...
            'xdir', 'reverse', ...
            'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
        myCM    = circshift(hsv, size(hsv, 1) / 2); % so that 0° is red
        colormap(ax1, myCM);
        axis off;
        
        % observations per angle
        ax2 = axes('units', 'centimeters', 'position', ax1.Position);
        maxN = max(NcircM);
        hold on;
        for iA = 1:size(dt.r.locDir.angleCenters, 1)
            tmpAngles   = transpose([0, dt.r.locDir.angleCenters(iA) - deg2rad(dt.r.locDir.angularRes / 2):...
                0.001:dt.r.locDir.angleCenters(iA) + deg2rad(dt.r.locDir.angularRes / 2), 0] + rotAngle);
            tmpRadii    = [0; repmat(NcircM(iA), numel(tmpAngles) - 2, 1); 0] ./ maxN;
            [x, y]      = pol2cart(tmpAngles, tmpRadii);
            patch(x, y, [0.7, 0.7, 0.7], ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(x, y, '-', ...
                'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        end        
        hold off;
        % adjust axes
        set(gca, ...
            'xdir', 'reverse', ...
            'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
        axis off;
    
        % save figure
        set(f, 'PaperPositionMode', 'auto');
        print(f, strcat(dt.r.paths.save, groups{iGroup}, '_', num2str(degSpace(iD)), 'DegSpace_NCircM_241119'), '-dtiff', '-r300');
    end
end

%% homogeneity of vector field for ELRPD cells vs. direction cells

% for each cell and each LRP, get the vector-field strength
vFS = nan(size(dt.r.allRes, 1), 1);
for iCell = 1:size(dt.r.allRes, 1)
    
    % mean preferred allocentric direction per LRP
    meanPDPerLRP    = nan(size(dt.r.allRes(iCell).dir_corrFR_perLRP, 1), 1);
    for iLRP = 1:size(dt.r.allRes(iCell).dir_corrFR_perLRP, 1)
        bNaN        = isnan(dt.r.allRes(iCell).dir_corrFR_perLRP(iLRP, :));
        if sum(~bNaN) >= 2
            meanPDPerLRP(iLRP)  = circ_mean(dt.direc.angleCenters(~bNaN), transpose(dt.r.allRes(iCell).dir_corrFR_perLRP(iLRP, ~bNaN)));
        end
    end
    
    % estimate vector field strength
    tmpMask          	= dt.r.locDir.LRPs2Exclude == false & ~isnan(meanPDPerLRP);
    if sum(tmpMask) > 0
        vFS(iCell, 1)   = circ_r(meanPDPerLRP(tmpMask));
    end
end

% check whether vector field strengths are higher for direction cells as
% compared to ELRPD cells
[~, p_vFS_dirVSELRPDCells, ~, stats_vFS_dirVSELRPDCells]    = ttest2(vFS(dt.bDirCell), vFS(dt.bELRPDCell));
fprintf('Vector-field strength for direction cells vs. ELRPD cells: %.3f vs. %.3f a.u. (t(%d) = %.3f, p = %.3f).\n', ...
    mean(vFS(dt.bDirCell)), mean(vFS(dt.bELRPDCell)), ...
    stats_vFS_dirVSELRPDCells.df, stats_vFS_dirVSELRPDCells.tstat, p_vFS_dirVSELRPDCells);

% mean, and standard error of the mean
m 	= [mean(vFS(dt.bDirCell)), mean(vFS(dt.bELRPDCell))];
ste = [std(vFS(dt.bDirCell)) ./ sqrt(sum(dt.bDirCell)), ...
    std(vFS(dt.bELRPDCell)) ./ sqrt(sum(dt.bELRPDCell))];

% create bar plot
f = figure('units', 'centimeters', 'position', [2, 5, 8, 8]);
axes('units', 'centimeters', 'position', [1.75, 2.25, 6, 5.5]);
hold on;
for iBar = 1:numel(m)
    b1 = bar(iBar, m(iBar), ...
        'FaceColor', [0, 0.5, 1]);
    if mod(iBar, 2) == 0
        set(b1, 'FaceColor', [0, 1, 0]);
    elseif mod(iBar, 3) == 0
        set(b1, 'FaceColor', [0.75, 0.75, 0.75]);
    end
    plot([iBar, iBar], [m(iBar) - ste(iBar), m(iBar) + ste(iBar)], ...
        'LineWidth', 2, 'Color', [0, 0, 0]);
end
% enhance axes
hold off;
set(gca, ...
    'xtick', 1:numel(m), 'xticklabel', '', ...
    'xlim', [0.3, numel(m) + 0.7], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
myAxes          = get(gca);
myXTickLabel    = {'Direction'; {'Egocentric', 'bearing'}};
for iLabel = 1:size(myXTickLabel, 1)
    text(myAxes.XTick(iLabel), -0.02, myXTickLabel{iLabel}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica', 'FontUnits', 'centimeters', 'FontSize', 0.5);
end
xl = xlabel('Cell type', ...
    'Position', [1.5, -0.22]);
yl = ylabel('Vector-field strength (a.u.)');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'VecFieldStrength_DirVSELRPDCells_20210409'), '-dtiff', '-r300');
    
%% overlap between ELRPD cells and place-like cells and direction cells

% overlap of different groups
A = [sum(dt.bDirCell), ...
    sum(dt.bPlaceCell), ...
    sum(dt.bELRPDCell)];
I = [sum(dt.bDirCell & dt.bPlaceCell), ...
    sum(dt.bDirCell & dt.bELRPDCell), ...
    sum(dt.bPlaceCell & dt.bELRPDCell), ...
    sum(dt.bDirCell & dt.bPlaceCell & dt.bELRPDCell)];

% create Venn diagram
f = figure('units', 'centimeters', 'position', [2, 5, 12, 8]);
venn(A, I, 'ErrMinMode', 'None', ...
    'FaceColor', {[0, 0.5, 1], [0, 0, 1], [0, 1, 0]}, ...
    'LineWidth', 2);
axis equal tight off;
set(gca, ...
    'ydir', 'reverse');
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'Venn_ELPRDCells_vs_PlaceCells_vs_DirCells_130919'), '-dtiff', '-r300');
