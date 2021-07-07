function [chi2stat, p] = myChiSquareTest(n11, n12, n21, n22)
%
% myChiSquareTest performs a Chi-squared test given the four entries
%   n11 n12
%   n21 n22
% from a contingency table.
%
% Use as: [chi2stat, p] = myChiSquareTest(n11, n12, n21, n22);
%
% See also: https://www.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
%
% Lukas Kunz, 2021

% basic populations
N1  = n11 + n21;
N2  = n12 + n22;

% pooled estimate of proportion
p0  = (n11 + n12) / (N1 + N2);

% expected counts under H0 (null hypothesis)
n110    = N1 * p0;
n120    = N2 * p0;

% chi-square test
observed    = [n11, N1 - n11, n12, N2 - n12];
expected    = [n110, N1 - n110, n120, N2 - n120];
chi2stat    = sum((observed - expected) .^ 2 ./ expected);
p           = 1 - chi2cdf(chi2stat, 1);
