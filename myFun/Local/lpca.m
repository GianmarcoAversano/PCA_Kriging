function [eigvec, n_eig, gamma, u_scores, eigenvalues, centroids] = lpca(nz_X_k, n_eigs, cent_crit, scal_crit, is_parallel)
% Centering and scaling criterion
if nargin < 3
    cent_crit = 1;
end
if nargin < 4
    scal_crit = 0;
end
if nargin < 5
    is_parallel = false;
end
% Number of clusters
k = length(nz_X_k);
n_vars = size(nz_X_k{1}, 2);
% Initialization of cell arrays
eigvec = cell(k, 1);
u_scores = cell(k, 1);
n_eig = cell(k, 1);
gamma = cell(k, 1);
eigenvalues = cell(k, 1);
centroids = cell(k, 1);
% Apply PCA in each cluster
if is_parallel
    parfor j = 1 : k
        % Center and scale, then do PCA
        [X, centroids{j}] = center(nz_X_k{j}, cent_crit);
        [X, gamma{j}] = scale(X, nz_X_k{j}, scal_crit);
        [modes, scores, eigenvalues{j}] = pca(X, 'Centered', false, 'Algorithm', 'svd'); 
        % Check n_eigs does not exceed the found number of modes
        n_modes = n_eigs;
        n_scores = n_eigs;
        if n_eigs > size(modes, 2)
            n_modes = size(modes, 2);
            n_scores = size(modes, 2);
        end
        % Outputs (and gamma)
        n_eig{j} = n_eigs;
        eigvec{j} = modes(:, 1:n_modes);
        u_scores{j} = scores(:, 1:n_scores);
    end
else
    for j = 1 : k
        % Center and scale, then do PCA
        [X, centroids{j}] = center(nz_X_k{j}, cent_crit);
        [X, gamma{j}] = scale(X, nz_X_k{j}, scal_crit);
        [modes, scores, eigenvalues{j}] = pca(X, 'Centered', false, 'Algorithm', 'svd'); 
        % Check n_eigs does not exceed the found number of modes
        n_modes = n_eigs;
        n_scores = n_eigs;
        if n_eigs > size(modes, 2)
            n_modes = size(modes, 2);
            n_scores = size(modes, 2);
        end
        % Outputs (and gamma)
        n_eig{j} = n_eigs;
        eigvec{j} = modes(:, 1:n_modes);
        u_scores{j} = scores(:, 1:n_scores);
    end
end
end



