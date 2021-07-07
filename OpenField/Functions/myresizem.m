function outMatrix = myresizem(inMatrix, augFac)
%
% MYRESIZEM resizes matrix "inMatrix" using the augmentation factor
% augFac. inMatrix can be an n x 1 vector or an n x m matrix.
%
% Output is the resized matrix "outMatrix".
%
% Lukas Kunz, 2021

% size of input matrix
[n, m]      = size(inMatrix);

% initialize output matrix
outMatrix   = cell(n, m);

% fill output matrix with values from 2D input matrix
if n > 1 && m > 1
    for i = 1:size(inMatrix, 1)
        for j = 1:size(inMatrix, 2)
            outMatrix{i, j}   = repmat(inMatrix(i, j), augFac, augFac);
        end
    end
end

% fill output matrix with values from 1D input matrix
if n > 1 && m == 1
    for i = 1:size(inMatrix, 1)
        for j = 1:size(inMatrix, 2)
            outMatrix{i, j}   = repmat(inMatrix(i, j), augFac, 1);
        end
    end
end

% unfold output matrix
outMatrix   = cell2mat(outMatrix);
