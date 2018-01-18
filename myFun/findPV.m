function [pv] = findPV(clusters, cen_crit, scal_crit, varargin)
%% Description

%% Input
% Optional arguments
if nargin < 2
    cen_crit = 1;
end
if nargin < 3
    scal_crit = 0;
end
% Check first argument
if ~iscell(clusters)
    clusters = {cluster};
end

%% Main
n_clust = length(clusters);
pv = cell(n_clust,1);
a_tol = 1e-16;
for ii = 1 : n_clust
    [~, ~, sig] = zscore(clusters{ii}, 0, 1);
    temp = center(clusters{ii}, cen_crit);
    temp = scale(temp, clusters{ii}, scal_crit);
    [modes, scores, eigens] = pca(temp, 'Centered', false);
    n_pc = size(modes, 2);
    y_dim = size(clusters{ii}, 2);
    loadings = modes * diag(eigens);
    for jj = 1 : y_dim
        loadings(jj,:) = loadings(jj,:) / (sig(jj) + a_tol);
    end
    pv{ii} = max(loadings, [], 1);
end
end
