function [clusters] = get_clusters(Y, idx, is_rows)
if nargin < 3
    is_rows = true;
end
n_clust = numel(unique(idx));
clusters = cell(n_clust, 1);
if is_rows
    for ii = 1 : n_clust
        clusters{ii} = Y(idx == ii, :);
    end
else
    for ii = 1 : n_clust
        clusters{ii} = Y(:, idx == ii);
    end
end
end
