
function [h,p,X2] = chi2cont(x,varargin)
% chi2cont chi-square test of contingency table
%   h = chi2cont(x) performs a chi-square test on the data in the
%   m-by-n contingency table x. The null hypothesis is that there is no difference
%   in the row variable distribution ('outcomes') between the columns 
%   ('treatments'). The result of the test is returned in h. h=1 indicates
%   a rejection of the null hypothesis at the 5% significance level.h=0
%   indicates that the null hypothesis can not be rejected at the 5%
%   significance level.
%
%   h = chi2cont(x,alpha) performs the test at the (100*alpha)%
%   significance level. The default when unspecified is alpha=0.05;
%   
%   [h,p] = chi2cont(...) returns the p value of the test. The p value is 
%   the probability, under the null hypothesis, of observing a value as 
%   extreme or more extreme of the chi-square test statistic.
%
%   [h,p,X2] = chi2cont(...) returns the chi-square test statistic.
%
%   Reference http://www.psychstat.missouristate.edu/introbook/sbk28m.htm
%
%   Mark Snaterse, January 22 2014
if isempty(varargin)
    alpha = 0.05;
elseif nargin>1
    error('Too many input arguments')
else
    alpha = varargin{1};
end
% Compute expectation and chi-square statistic, and determine p value.
e = sum(x,2)*sum(x)/sum(x(:));
X2 = (x-e).^2./e;
X2 = sum(X2(:));
df = prod(size(x)-[1 1]);
p = 1-chi2cdf(X2,df);
h = double(p<=alpha);
