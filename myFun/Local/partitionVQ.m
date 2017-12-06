function [nz_X_k, nz_idx_clust, k, varargout] = partitionVQ(X, idx, varargin)
[n_samples, n_variables] = size(X); % n_samples; n_variables
threshold = n_variables;
% Remove cluster if n_samples < threshold
isremove = true; % Default
if ~isempty(varargin)
    isremove = varargin{1}; % User-provided input
end
if ~isa(isremove, 'logical')
    isremove = false; % User was trying to provide an input: interpret as false
end
k = numel(unique(idx)); % Number of clusters = number of unique integers in IDX
idx_clust = cell(k, 1); % Initialize
n_points = zeros(k, 1); % Initialize
for j = 1 : k
    idx_clust{j} = find(idx == j); % Rows of X belonging to cluster j
    n_points(j) = size(idx_clust{j}, 1); % Clusters' population
    if isremove
        % Remove clusters that do not have enough population
        if (n_points(j) < threshold)
            fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
        end
    end
end
if ~isremove
    nz_idx = find(n_points > -1); % Every cluster has a population > -1: no clusters will be removed
else
    nz_idx = find(n_points > threshold); % Only this clusters will be kept
end
k = size(nz_idx, 1); % New number of clusters
nz_X_k = cell(k, 1); % Initialize
nz_idx_clust = cell(k, 1); % Initialize 
for j = 1 : k
    nz_X_k{j} = zeros(n_points(j), n_variables); % Cluster j
    nz_idx_clust{j} = idx_clust{nz_idx(j)}; % New IDX and IDX_CLUST: in case the number of clusters changed 
    nz_X_k{j} = X(nz_idx_clust{j}, :); % Cluster j
end
% Optional output
if nargout > 3
    idx = removeHoles(nz_idx, true);
    varargout{1} = idx;
end
end



 
       