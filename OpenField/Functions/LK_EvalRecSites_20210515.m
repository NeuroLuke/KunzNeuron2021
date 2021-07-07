function [] = LK_EvalRecSites_20210515(dt)
%
% LK_EvalRecSites_20210515 evaluates the influence of recording sites.
%
% Input is structure with fields
%   r        	--> multiple results from the ELRPD-cell analysis
%   bELRPDCell 	--> boolean vector to indicate ELRPD cells
%
% Lukas Kunz, 2021

%% epileptic region

% report
fprintf('\n----- Evaluation of ELRPD cells when excluding epileptic brain regions.\n');

% epileptic regions
epilepticReg    = cell(size(dt.r.allRes, 1), 1);
for iCell = 1:size(dt.r.allRes, 1)
    
    % load subject information
    s   = load(strcat(dt.r.paths.info, dt.r.subjects{dt.r.allRes(iCell).idx(1)}, '\subjectdata_180619.mat'));    
    
    % for this wire, extract information about epileptic regions
    logIdx             	= any(cell2mat(s.subjectdata.micro2macro(:, 1)) == dt.r.allRes(iCell).idx(2), 2);
    epilepticReg{iCell} = s.subjectdata.micro2macro{logIdx, 6};
end

% exclude wires involved in the seizure origin
bNonIctal  = cellfun(@isempty, regexp(epilepticReg, 'ictal'));
fprintf('Total number of cells not recorded from "ictal" wires: %d.\n', sum(bNonIctal));

% report number of cells when excluding ictal regions
fprintf('Number of ELRPD cells when excluding "ictal" wires: %d, %.3f%%, binomial P = %.3f.\n', ...
    sum(dt.bELRPDCell & bNonIctal), ...
    100 * sum(dt.bELRPDCell & bNonIctal) / sum(bNonIctal), ...
    myBinomTest(sum(dt.bELRPDCell & bNonIctal), sum(bNonIctal), 0.05));

%% hemisphere

% differentiate between right and left hemisphere
allMNI              = cell2mat({dt.r.allRes.wireMNI}');
bRightHemisphere    = allMNI(:, 1) > 0;

% contingency table
n   = [sum(dt.bELRPDCell & bRightHemisphere), sum(dt.bELRPDCell & ~bRightHemisphere); ...
    sum(~dt.bELRPDCell & bRightHemisphere), sum(~dt.bELRPDCell & ~bRightHemisphere)];
[X, p] = myChiSquareTest(n(1, 1), n(1, 2), n(2, 1), n(2, 2));
fprintf('Are there more ELRPD cells in one hemisphere? X2 = %.3f, p = %.3f.\n', X, p);

% report
fprintf('Number of ELRPD cells in the right hemisphere: %d (of %d, %.3f%%).\n', sum(dt.bELRPDCell & bRightHemisphere), sum(bRightHemisphere), ...
    sum(dt.bELRPDCell & bRightHemisphere) / sum(bRightHemisphere) * 100);
fprintf('Number of ELRPD cells in the left hemisphere: %d (of %d, %.3f%%).\n', sum(dt.bELRPDCell & ~bRightHemisphere), sum(~bRightHemisphere), ...
    sum(dt.bELRPDCell & ~bRightHemisphere) / sum(~bRightHemisphere) * 100);
