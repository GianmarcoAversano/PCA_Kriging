function X = uncenter(X, mean_vals)
%% Description:
% X: (n_obs x n_vars)
% mean_vals: array of mean values
%

%% Main
[rows, ~] = size(X);
X = X + repmat(mean_vals(:)', rows, 1);

end


