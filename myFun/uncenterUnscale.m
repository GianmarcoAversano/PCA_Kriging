function [varargout] = uncenterUnscale(X, mean_column, scaling_factors, varargin)
%% Description
% X - data, n_vars x n_obs
% 

%% Input
% Check size of input
[r, c] = size(mean_column);
if r > 1 && c > 1
    error('Column of mean values has wrong size.');
end
[r, c] = size(scaling_factors);
if r > 1 && c > 1
    error('Column of scaling factors has wrong size.');
end
% Reshape input
if length(mean_column) == 1
    mean_column = zeros(size(X,1),1) + mean_column;
end
if length(scaling_factors) == 1
    scaling_factors = zeros(size(X,1),1) + scaling_factors;
end

%% Main
% Relavant dimension
[n_rows, n_cols] = size(X);
lm = length(mean_column);
ls = length(scaling_factors);
% Check sizes
if (lm ~= n_rows) || (ls ~= n_rows)
    error('[uncenterUnscale:] size(X,1) = %i; length(m) = %i; length(sc) = %i', n_rows, lm, ls);
end
% Diagonal matrix of scaling factors
D = spdiags(scaling_factors(:), 0, n_rows, n_rows);
% Unscale
X = D * X;
% Uncenter
M = repmat(mean_column(:), 1, n_cols);
% Decoded data
X = M + X;

%% Output
if nargout > 0
    varargout{1} = X;
end

end