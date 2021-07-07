function [] = LK_TH_EvalELRPDCells_20210515(dt)
%
% LK_TH_EvalELRPDCells_20210515 evaluates ELRPD cells in Treasure Hunt.
%
% Input is a structure with multiple fields related to the ELRPD cell
% analysis.
% 
% Lukas Kunz, 2021

% report
fprintf('\n----- Evaluating ELRPD cells.\n');
rng(1);

% font size
myFontSize  = 0.5;

% rotate egocentric directions so that "ahead" is pointing "north"
rotAngle    = deg2rad(90);

%% number of ELRPD cells per patient

% report
fprintf('\nEvaluating the number of ELRPD cells per patient and session.\n');

% unit indices
allUnitIdx              = cell2mat({dt.r.allRes.idx}');

% number of ELRPD cells per subject and session
numCellsPerSubj         = cell(size(dt.r.subjects, 1), 1);
numELRPDCellsPerSubj    = cell(size(dt.r.subjects, 1), 1);
for iSubj = 1:size(dt.r.subjects, 1)
    
    % subject's sessions
    sessIdx             = unique(allUnitIdx(allUnitIdx(:, 1) == iSubj, 2));
    
    % number of cells per session
    tmpNumCells         = nan(size(sessIdx, 1), 1);
    tmpNumELRPDCells    = nan(size(sessIdx, 1), 1);
    for iSess = 1:numel(sessIdx)
        tmpNumCells(iSess, 1)       = sum(allUnitIdx(:, 1) == iSubj & allUnitIdx(:, 2) == sessIdx(iSess));
        tmpNumELRPDCells(iSess, 1)  = sum(allUnitIdx(:, 1) == iSubj & allUnitIdx(:, 2) == sessIdx(iSess) & dt.bELRPDCell == true);
    end
    numCellsPerSubj{iSubj, 1}       = tmpNumCells;
    numELRPDCellsPerSubj{iSubj, 1}  = tmpNumELRPDCells;
end

% report
fprintf('Number of sessions with at least one ELRPD cell: %d (out of %d).\n', sum(cell2mat(numELRPDCellsPerSubj) > 0), numel(cell2mat(numELRPDCellsPerSubj)));
fprintf('Average number of ELRPD cells per session: %.3f +/- %.3f.\n', ...
    mean(cell2mat(numELRPDCellsPerSubj)), std(cell2mat(numELRPDCellsPerSubj)) / sqrt(numel(cell2mat(numELRPDCellsPerSubj))));

% number of unique subjects with at least one ELRPD cell
fprintf('Number of unique subjects with at least one ELRPD cell: %d (out of %d).\n', sum(cellfun(@(x) sum(x), numELRPDCellsPerSubj) > 0), size(numELRPDCellsPerSubj, 1));

%% evaluate distance of COMs (= reference points) to environmental center

% report
fprintf('\nEvaluating the center distance of reference points (i.e., reference-field COMs).\n');

% COMs
allCOMxz                = cell2mat({dt.r.allRes.locdir_COMxz}');
allCOMD2Ctr             = pdist2(allCOMxz, dt.r.param.arenaCtr); % center distances of COMs
ELRPDCells_allCOMD2Ctr  = allCOMD2Ctr(dt.bELRPDCell, 1); % center distances of COMs, ELRPD cells only

% step through each center-distance bin and estimate significance
rTotal          = dt.r.param.arenaRadius; % radius of entire arena
areaTotal       = pi * rTotal ^ 2; % area of entire arena
rEdges          = 0:2.5:rTotal;
rCenters        = movmean(rEdges, 2, 'endpoints', 'discard');
numCOMsPerR     = nan(numel(rCenters), 1); % number of COMs per distance
numCOMsPerRExp  = nan(numel(rCenters), 1);
pCOMsPerR       = nan(numel(rCenters), 1); % p-values
for iR = 1:numel(rCenters)
    
    % observed number of reference points
    innerR                  = rEdges(iR);
    outerR                  = rEdges(iR + 1);
    numCOMsPerR(iR, 1)      = sum(ELRPDCells_allCOMD2Ctr >= innerR & ELRPDCells_allCOMD2Ctr < outerR);
    
    % expected number of reference points
    numCOMsPerRExp(iR, 1)   = numel(ELRPDCells_allCOMD2Ctr) * (pi * outerR ^ 2 - pi * innerR ^ 2) / areaTotal;
    
    % binomial test
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
fprintf('The border between the two groups (n = %d and n = %d) is at %.3f.\n', ...
    sum(ELRPDCells_allCOMD2Ctr_Group == 1), sum(ELRPDCells_allCOMD2Ctr_Group == 2), borderCtrPeri);

%% figure showing the COM distances to the environmental center

% histogram for distance of COMs to the environmental center
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
axes('units', 'centimeters', 'position', [1.55, 1.75, 6, 6]);
hold on;
% histograms of both groups
H1 = histogram(ELRPDCells_allCOMD2Ctr(ELRPDCells_allCOMD2Ctr_Group == 1), rEdges, ...
    'facecolor', rgb('darkgreen'), 'edgecolor', 'none', 'facealpha', 1); % center reference points
H2 = histogram(ELRPDCells_allCOMD2Ctr(ELRPDCells_allCOMD2Ctr_Group == 2), rEdges, ...
    'facecolor', rgb('limegreen'), 'edgecolor', 'none', 'facealpha', 1); % periphery reference points
H3 = stairs(rEdges, [numCOMsPerRExp; numCOMsPerRExp(end)], ...
    'Color', [0.5, 0.5, 0.5]);
axis square;
set(gca, ...
    'xlim', [min(rEdges) - 1, max(rEdges) + 1], 'xtick', [0, 25, 50], ...
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
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'ELRPDCells_COMxzDistance2Center_120919'), '-dtiff', '-r300');

%% evaluate COMs

% report
fprintf('\nEvaluating the spatial distribution of reference points (i.e., reference-field COMs).\n');

% reference points of ELRPD cells
ELRPDCells_allCOMxz = allCOMxz(dt.bELRPDCell, :);

% plot centers of mass within the environment
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
axes('units', 'centimeters', 'position', [1.15, 1.15, 6, 6]);
hold on;
% boundary
plot(dt.r.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* dt.r.param.arenaRadius, dt.r.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* dt.r.param.arenaRadius, '-', ...
    'Color', [0, 0, 0], 'LineWidth', 5);
% center/periphery boundary
plot(dt.r.param.arenaCtr(1) + cos(0:0.001:(2*pi)) .* borderCtrPeri, dt.r.param.arenaCtr(2) + sin(0:0.001:(2*pi)) .* borderCtrPeri, ':', ...
    'Color', [0, 0, 0], 'LineWidth', 2);
% reference points in the center and periphery
scatter(ELRPDCells_allCOMxz(ELRPDCells_allCOMD2Ctr_Group == 1, 1), ...
    ELRPDCells_allCOMxz(ELRPDCells_allCOMD2Ctr_Group == 1, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('darkgreen'), 'MarkerFaceColor', rgb('darkgreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1);
scatter(ELRPDCells_allCOMxz(ELRPDCells_allCOMD2Ctr_Group == 2, 1), ...
    ELRPDCells_allCOMxz(ELRPDCells_allCOMD2Ctr_Group == 2, 2), 50, 'o', ...
    'MarkerEdgeColor', rgb('limegreen'), 'MarkerFaceColor', rgb('limegreen'), 'MarkerFaceAlpha', 0.4, 'LineWidth', 1);
hold off;
set(gca, ...
    'xdir', 'reverse', 'ydir', 'reverse', ... % flip x- and y-axis
    'xlim', [min(dt.r.locDir.xEdges), max(dt.r.locDir.xEdges)], 'xtick', [min(dt.r.locDir.xCenters), max(dt.r.locDir.xCenters)], ...
    'ylim', [min(dt.r.locDir.zEdges), max(dt.r.locDir.zEdges)], 'ytick', [min(dt.r.locDir.zCenters), max(dt.r.locDir.zCenters)], ...
    'tickdir', 'out', 'ticklength', [0, 0], ...
    'box', 'on');
xl = xlabel('x (vu)', ...
    'units', 'normalized', 'position', [0.5, -0.025, 0]);
yl = ylabel('y (vu)', ...
    'units', 'normalized', 'position', [-0.025, 0.5, 0]);
ytickangle(90);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'ELRPDCells_allCOMxz_120919'), '-dtiff', '-r300');

%% distribution of preferred ELRPDs

% for each cell, get preferred ELRPD
prefELRPD   = nan(size(dt.r.allRes, 1), 1);
for iCell = 1:size(prefELRPD, 1)
    % if the cell has a COM, extract preferred ELRPD
    if ~isempty(dt.r.allRes(iCell).locdir_COMxzClosestLRP)
        thisCOMLRP          = all(dt.r.locDir.LRPs == dt.r.allRes(iCell).locdir_COMxzClosestLRP, 2);
        prefELRPD(iCell, 1) = circ_mean(transpose(dt.r.locDir.angleCenters), ...
            transpose(dt.r.allRes(iCell).all_locdir_corrFR(thisCOMLRP == true, :)));
    end
end

% examine different groups of cells and different angular spaces
groups      = {'ELRPDCellsCtrCOM'; 'ELRPDCellsPeriCOM'};
degSpace    = [360; 180]; % different angular spaces

% loop through groups
all_Rz_diffGroupsDiffSpaces = nan(size(groups, 1), size(degSpace, 1));
all_Rp_diffGroupsDiffSpaces = nan(size(groups, 1), size(degSpace, 1));
for iGroup = 1:size(groups, 1)
    
    % loop through different angular spaces
    for iD = 1:size(degSpace, 1)
        
        % group of cells to analyze
        if strcmp(groups{iGroup}, 'ELRPDCellsCtrCOM')
            bCellMask   = dt.bELRPDCell & allCOMD2Ctr < borderCtrPeri;
        elseif strcmp(groups{iGroup}, 'ELRPDCellsPeriCOM')
            bCellMask   = dt.bELRPDCell & allCOMD2Ctr > borderCtrPeri;
        end
        
        % Rayleigh test in specific angular space
        [pval, z] = circ_rtest(prefELRPD(bCellMask) .* (360 / degSpace(iD)));
        fprintf('Rayleigh test for group "%s" in %d°-space: z = %.3f, P = %.3f.\n', ...
            groups{iGroup}, degSpace(iD), z, pval);
        all_Rz_diffGroupsDiffSpaces(iGroup, iD) = z;
        all_Rp_diffGroupsDiffSpaces(iGroup, iD) = pval;
        
        % number of observations per angular bin
        tmpF    = figure('visible', 'off');
        myP     = polarhistogram(prefELRPD(bCellMask) .* (360 / degSpace(iD)), dt.r.locDir.angleEdges);
        NcircM  = myP.Values;
        close(tmpF);
        
        %--- create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
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
        tmpx = (-pi:0.001:pi) + rotAngle; % rotate so that red is pointing "north"
        tmpc = linspace(0, 1, length(tmpx));
        patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
            [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
            'EdgeColor', 'interp');
        % plot angle text
        if iGroup == 1
            t1 = text(1.125, 0, 'R');
            t2 = text(0, 1.15, 'A');
            t3 = text(-1.125, 0, 'L');
            t4 = text(0, -1.15, 'B');
            set([t1, t2, t3, t4], ...
                'FontUnits', 'centimeters', 'FontSize', myFontSize, ...
                'HorizontalAlignment', 'center');
        end
        % indicate maximum
        text(0.175, 0.15, ['max ', num2str(max(NcircM), '%.0f')], ...
            'units', 'normalized', ...
            'horizontalalignment', 'right', 'verticalalignment', 'top', ...
            'fontunits', 'centimeters', 'fontsize', myFontSize);
        hold off;
        % adjust x-axis
        set(gca, ...
            'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
        myCM    = flipud(circshift(hsv, size(hsv, 1) / 2));
        colormap(ax1, myCM);
        axis off;
        
        % observations per angle
        ax2 = axes('units', 'centimeters', 'position', ax1.Position);
        maxN = max(NcircM);
        hold on;
        for iA = 1:numel(dt.r.locDir.angleCenters)
            % plot angles and corresponding data
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
        set(ax2, ...
            'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
        axis off;
    
        % save figure
        set(f, 'PaperPositionMode', 'auto');
        print(f, strcat(dt.r.paths.save, groups{iGroup}, '_', num2str(degSpace(iD)), 'DegSpace_NCircM_241119'), '-dtiff', '-r300');
    end
end
