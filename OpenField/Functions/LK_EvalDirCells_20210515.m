function [] = LK_EvalDirCells_20210515(dt)
%
% LK_EvalDirCells_20210515 evaluates direction cells.
%
% Lukas Kunz, 2021

fprintf('\nDirection-cell evaluation.\n');

%% clustering of preferred directions

% for each cell, determine circular mean of the firing rate map, i.e., the
% preferred direction
circM   = nan(size(dt.r.allRes, 1), 1);
for iCell = 1:size(dt.r.allRes, 1)
    circM(iCell, 1) = circ_mean(dt.r.direc.angleCenters, dt.r.allRes(iCell).dir_corrFR);
end

% evaluate circular uniformity of preferred directions of direction cells
[pval_Rayleigh_circM, z_Rayleigh_circM] = circ_rtest(circM(dt.bDirCell));
fprintf('Non-uniformity of preferred directions of direction cells? z = %.3f, P = %.3f.\n', ...
    z_Rayleigh_circM, pval_Rayleigh_circM);

