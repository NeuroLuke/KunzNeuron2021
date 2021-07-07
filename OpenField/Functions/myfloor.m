function outNum = myfloor(inNum, precision)
%
% MYFLOOR is identical to floor with the additional option to specify a
% precision to floor towards a specific decimal place. The precision is
% specified as a power of 10 (with positive precision values leading to
% more precise flooring).
%
% Lukas Kunz, 2021

prec    = 10 .^ precision;
outNum  = floor(inNum .* prec) ./ prec;