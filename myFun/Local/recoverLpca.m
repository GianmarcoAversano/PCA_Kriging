function lpca_X = recoverLpca(local_modes, local_scores, centroids, idx, mean_val, sig_val)
% Dummy check
if length(local_modes) ~= length(local_scores)
    error('No coherence found.');
end
% Get info
n_clust = length(local_modes);
n_vars = size(local_modes{1}, 1);
n_samples = 0;
for ii = 1 : length(local_scores)
    n_samples = n_samples + size(local_scores{ii}, 1);
end
% Initialize matrix
lpca_X = zeros(n_samples, n_vars);
% Recover data from Local PCA
for ii = 1 : n_clust
    % Recover data from one cluster
    I = (idx == ii);
    lpca_X(I,:) = local_scores{ii} * local_modes{ii}' + ...
        repmat(centroids(ii,:), size(lpca_X(I,:),1), 1);
end
% This data is still centered-scaled. We have recovered the centered-scaled
% data that was passed as input to the clustering/local PCA procedure.
if nargin > 4
    lpca_X = uncenterUnscale(lpca_X, mean_val, sig_val);
end
end

