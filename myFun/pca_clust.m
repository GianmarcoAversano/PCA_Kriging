function [idx, clusters, clusters_idx, T] = pca_clust(Y, n_pc, varargin)
if nargin < 2
    n_pc = min(size(Y)); % Number of PCs used for the splitting
end
y_dim = size(Y,2);
[modes, scores, eigens] = pca(Y, 'Centered', false, 'Algorithm', 'svd');
n_scores = size(scores,2);
if n_pc > n_scores
    n_pc = n_scores;
end
M = size(scores, 1); % Number of observations
% Evaluate the affinities of each variable
T = zeros(n_pc, y_dim);
for ii = 1 : n_pc
    T(ii,:) = modes(:,ii)';
end
% Assing each variable to a cluster: assignment depends on the affinities
% to the PCs
[~, idx] = max(T, [], 1);
idx = removeHoles(idx, true);
n_clust = max(idx); % Actual number of clusters found
clusters = cell(n_clust,1);
clusters_idx = cell(n_clust,1);
for ii = 1 : n_clust
    clusters{ii} = Y(:, idx == ii);
    clusters_idx{ii} = idx(idx == ii);
end
end
