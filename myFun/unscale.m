function X = unscale(X, gamma)
%% Description:
% X: (n_obs x n_vars)
% gamma: array of scaling factors
%

%% Main
[~, columns] = size(X);
if length(gamma) ~= columns
    error('Array GAMMA must be as long as columns of X.');
end
parfor ii = 1 : columns
    X(:,ii) = X(:,ii) * gamma(ii);
end

end


