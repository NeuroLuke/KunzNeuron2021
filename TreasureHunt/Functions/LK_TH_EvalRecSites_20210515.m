function [] = LK_TH_EvalRecSites_20210515(dt)
%
% LK_TH_EvalRecSites_20210515 evaluates the recording sites.
%
% Input is a structure with multiple fields.
%
% Lukas Kunz, 2021

% report
fprintf('\n----- Evaluation of recording sites.\n');

%% number of cells as a function of hemisphere

% epileptic regions and hemisphere
allLabels       = cell(size(dt.r.allRes, 1), 1);
epilepticReg    = cell(size(dt.r.allRes, 1), 1);
allHemi         = cell(size(dt.r.allRes, 1), 1);
for iCell = 1:size(dt.r.allRes, 1)
    
    % load subject information
    s                   = load(strcat(dt.r.paths.info, dt.r.subjects{dt.r.allRes(iCell).idx(1)}, '\subjectdata.mat'));    
    
    % for this wire, extract information about epileptic region
    logIdx              = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == dt.r.allRes(iCell).idx(3), 2);
    allLabels{iCell}    = s.subjectdata.micro2macro{logIdx, 2};
    allHemi{iCell}      = s.subjectdata.micro2macro{logIdx, 4};
    epilepticReg{iCell} = s.subjectdata.micro2macro{logIdx, 6};
end

% cells in right vs. left hemisphere
bRightHemi  = strcmp(allHemi, 'right');
fprintf('\nTotal number of cells: %d.\n', size(bRightHemi, 1));
fprintf('Number of cells recorded from the right hemisphere: %d.\n', sum(bRightHemi));
fprintf('Number of cells recorded from the left hemisphere: %d.\n', sum(~bRightHemi));

% non-ELRPD cells in the right vs. left hemisphere
fprintf('Number of non-ELRPD cells recorded from the right hemisphere: %d.\n', sum(~dt.bELRPDCell & bRightHemi));
fprintf('Number of non-ELRPD cells recorded from the left hemisphere: %d.\n', sum(~dt.bELRPDCell & ~bRightHemi));

% ELRPD cells in right vs. left hemisphere
fprintf('Number of ELRPD cells recorded from the right hemisphere: %d (%.3f%% of %d cells).\n', sum(dt.bELRPDCell & bRightHemi), ...
    100 * sum(dt.bELRPDCell & bRightHemi) / sum(bRightHemi), sum(bRightHemi));
fprintf('Number of ELRPD cells recorded from the left hemisphere: %d (%.3f%% of %d cells).\n', sum(dt.bELRPDCell & ~bRightHemi), ...
    100 * sum(dt.bELRPDCell & ~bRightHemi) / sum(~bRightHemi), sum(~bRightHemi));

%% evaluation whether more ELRPD cells are recorded from one hemisphere

% chi-squared test
n = [sum(dt.bELRPDCell & bRightHemi), sum(dt.bELRPDCell & ~bRightHemi); ...
    sum(~dt.bELRPDCell & bRightHemi), sum(~dt.bELRPDCell & ~bRightHemi)];
[X, p] = myChiSquareTest(n(1, 1), n(1, 2), n(2, 1), n(2, 2));
fprintf('Overlap between ELRPD cells (yes vs. no) and hemisphere (right vs. left)? X2 = %.3f, p = %.3f.\n', X, p);
disp(n);

%% number of cells as a function of epileptic vs. non-epileptic region

% report
fprintf('\nEvaluation of ELRPD cells when excluding epileptic brain regions.\n');

% exclude wires involved in the seizure origin
bNonIctal  = cellfun(@isempty, regexp(epilepticReg, 'ictal'));
fprintf('Total number of cells *not* recorded from "ictal" wires: %d.\n', sum(bNonIctal));

% report number of cells when excluding ictal regions
fprintf('Number of ELRPD cells when excluding "ictal" wires: %d, %.3f%%, binomial P = %.3f.\n', ...
    sum(dt.bELRPDCell & bNonIctal), ...
    100 * sum(dt.bELRPDCell & bNonIctal) / sum(bNonIctal), ...
    myBinomTest(sum(dt.bELRPDCell & bNonIctal), sum(bNonIctal), 0.05));

