function LK_EvalObjCells_090120(dt)
%
% LK_EvalObjCells_090120 evaluates object cells.
%
% Lukas Kunz, 2021

%% overlap between object cells and ELRPD cells

% overlap between object cells and ELRPD cells
A   = [sum(dt.bObjCell), sum(dt.bELRPDCell)];
I   = sum(dt.bObjCell & dt.bELRPDCell);

% venn diagram
f = figure('units', 'centimeters', 'position', [2, 2, 12, 8]);
venn(A, I, ...
    'FaceColor', {rgb('orange'), [0.5, 1, 0]}, 'LineWidth', 2);
axis equal tight off;
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'Venn_ObjCells_vs_ELPRDCells_130919'), '-dtiff', '-r300');

%% number of preferred objects

% identify preferred objects
sigBins     = cell2mat({dt.r.allRes.objBins})';

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 4, 4]);
h1 = histogram(sum(sigBins(dt.bObjCell, :), 2), ...
    'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', [1, 1, 1]);
tmpAx = get(gca);
set(gca, ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'xtick', movmean(h1.BinEdges, 2, 'endpoints', 'discard'), ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('# pref. objects');
yl = ylabel('Cell count');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.4);
% save image
set(gcf, 'PaperPositionMode', 'auto');
print(f, '-dtiff', ...
    strcat(dt.r.paths.save, 'ObjCells_numSigObj_20200812.tiff'), '-r300');

%% examine whether multiple preferred objects arise due to a spatial relationship of their locations

% reset rng
rng(444);

% preallocate
thisNumSurrogates   = 401; % number of surrogates
D_prefObjLocs       = nan(size(dt.r.allRes, 1), 1); % empirical distances
D_prefObjLocsSurro  = nan(size(dt.r.allRes, 1), thisNumSurrogates); % surrogate distances

% loop through cells
for iCell = 1:size(dt.r.allRes, 1)
    
    % object locations
    objLocs         = unique(dt.r.allBeh(dt.r.allRes(iCell).idx(1)).trials(:, [2, 9, 10]), 'stable', 'rows');
    objLocs         = objLocs(all(~isnan(objLocs), 2), :);
    [~, I]          = sort(objLocs(:, 1));
    objLocs         = objLocs(I, 2:3); % objects in correct order (0:7)
    
    % distance between preferred objects
    prefObjLocs     = objLocs(dt.r.allRes(iCell).objBins, :);
    if ~isempty(prefObjLocs)
        D_prefObjLocs(iCell, 1) = mean(pdist(prefObjLocs));
        
        % surrogate distances
        for iSurro = 1:thisNumSurrogates
            
            % randomly select n object locations
            randIdx                             = datasample(1:size(objLocs, 1), sum(dt.r.allRes(iCell).objBins), 'replace', false);
            D_prefObjLocsSurro(iCell, iSurro)   = mean(pdist(objLocs(randIdx, :)));
        end
    end
end

% report outcome
fprintf('Are the locations of preferred objects closer together in space than random objects?\n');
tmpRank = sum(nanmean(D_prefObjLocs(dt.bObjCell)) < nanmean(D_prefObjLocsSurro(dt.bObjCell, :))) / thisNumSurrogates;
fprintf('Rank of empirical value within (i.e., smaller than) surrogates = %.3f (p = %.3f).\n', ...
    tmpRank, 1 - tmpRank);
fprintf('Number of object cells contributing to this analysis: %d.\n', sum(~isnan(D_prefObjLocs(dt.bObjCell))));

%% create figure for statistical evaluation of spatial proximity of preferred objects

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7, 8]);
axes('units', 'centimeters', 'position', [1, 2.24682910989324, 5.25, 5.14430630677342]);
hold on;
myH = histogram(nanmean(D_prefObjLocsSurro(dt.bObjCell, :)), 2500:100:4500, ...
    'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
plot([nanmean(D_prefObjLocs(dt.bObjCell, :)), nanmean(D_prefObjLocs(dt.bObjCell, :))], [0, max(myH.Values)], ...
    'Color', [1, 0, 0], 'LineWidth', 2);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [2500, 4500], 'xtick', [2500, 4500], ...
    'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel({'Distance between', 'pref. objects (vu)'});
yl = ylabel('Count', ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(f, 'PaperPositionMode', 'auto');
print(f, strcat(dt.r.paths.save, 'ObjCells_DistBetwPrefObj_080320'), '-dtiff', '-r300');
