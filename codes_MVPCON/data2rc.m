function [rc_w,rc_b] = data2rc(Ya,Yb,metric,corrType);
% This function computes the representational connectivity from the
% multivariate time series of two regions. 
% It computes both the wihtin- and between-run RDM correlations for the two
% regions.
% Hamed Nili

if nargin < 3
    metric = 'correlation';
end
if ~exist('corrType','var')
    corrType = 'Pearson';
end

nruns = numel(Ya); % number of runs

for r=1:nruns
    rdms1(:,r)    = pdist(Ya{r},metric)';
    rdms2(:,r)    = pdist(Yb{r},metric)';
end
corrmat  = corr(rdms1,rdms2,'type',corrType);
rc_w = mean(diag(corrmat));
rc_b = mean(corrmat(logical(1-eye(nruns))));
