function [idx, clusters, clusters_idx] = hclust_pc(Y, ktot, varargin)
%% Description

%% Input
% Total number of clusters to find
if nargin < 2
    ktot = 2; 
end

%% Main
y_dim = size(Y, 2);
% Initialization
n_clust = 1;
clusters = {Y};
idx = ones(y_dim, 1);
% Split each cluster until max num of clusters is reached
cond = true;
while(cond)
    % Split each cluster
    idx_c = cell(n_clust, 1);
    k = zeros(n_clust,1);
    for ii = 1 : n_clust
        [idx_c{ii}, ~, c] = pca_clust(clusters{ii});
        k(ii) = length(c);
    end
    % Re-set global idx
    offset = 0;
    idx_old = idx;
    for ii = 1 : n_clust
        I = (idx_old == ii);
        idx(I) = idx_c{ii} + offset;
        % Update offset
        offset = max(idx);
    end
    idx = removeHoles(idx, true);
    % Update clusters
    n_clust_old = n_clust;
    n_clust = max(idx); % Update number of clusters found
    clusters = cell(n_clust,1);
    clusters_idx = cell(n_clust,1);
    for ii = 1 : n_clust
        clusters{ii} = Y(:, idx == ii);
        clusters_idx{ii} = idx(idx == ii);
    end
    % Convergence checks
    if n_clust <= n_clust_old
        % Check the number of clusters has increased
        cond = false;
    else
        % Check goal has been reached
        cond = (n_clust < ktot);
    end
end
end





%% OLD
function [idx, clusters, clusters_idx] = hclust_pc_old(Y, ktot, varargin)
if nargin < 2
    ktot = 2; % Total number of clusters to find
end
y_dim = size(Y,2);
% First split
[idx, clusters, clusters_idx] = pca_split(Y);
n_clust = length(clusters);
% Split each cluster until max num of clusters is reached
cond = true;
while(cond)
    % Split each cluster
    idx_c = cell(n_clust, 1);
    k = zeros(n_clust,1);
    for ii = 1 : n_clust
        [idx_c{ii}, ~, c] = pca_split(clusters{ii});
        k(ii) = length(c);
    end
    % Re-set global properties
    offset = 0;
    idx_old = idx;
    for ii = 1 : n_clust
        I = (idx_old == ii);
        idx(I) = idx_c{ii} + offset;
        % Update offset
        offset = max(idx);
    end
    idx = removeHoles(idx, true);
    n_clust_old = n_clust;
    n_clust = max(idx); % Update number of clusters found
    clusters = cell(n_clust,1);
    clusters_idx = cell(n_clust,1);
    for ii = 1 : n_clust
        clusters{ii} = Y(:, idx == ii);
        clusters_idx{ii} = idx(idx == ii);
    end
    % Convergence checks
    if n_clust <= n_clust_old
        % Check the number of clusters has increased
        cond = false;
    else
        % Check goal has been reached
        cond = (n_clust < ktot);
    end
end
end




