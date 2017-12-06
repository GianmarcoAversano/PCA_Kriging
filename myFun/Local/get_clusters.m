function [clusters] = get_clusters(Y, idx)
n_clust = numel(unique(idx));
clusters = cell(n_clust, 1);
for ii = 1 : n_clust
    clusters{ii} = Y(idx == ii, :);
end
end
