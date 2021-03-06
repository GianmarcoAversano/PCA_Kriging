function z = center_scale(X, mu, sigma, varargin)
% [] is a special case for std and mean, just handle it out here.
% Preference is: X (obs x var)
if isempty(X)
    z = []; 
    return
end
% Check mu and sigma are vectors
if size(mu,1) ~= 1 && size(mu,2) ~= 1
    error('MU must be a vector.');
end
if size(sigma,1) ~= 1 && size(sigma,2) ~= 1
    error('SIGMA must be a vector.');
end
% Tolerance
a_tol = 1e-16;
if ~isempty(varargin)
    a_tol = varargin{1};
end
% Get sizes
[rows, cols] = size(X);
% Subtract the mean
if length(mu) == rows
    z = X - repmat(mu(:), 1, cols);
elseif length(mu) == cols
    z = X - repmat(mu(:)', rows, 1);
end
% Scale
if length(sigma) == rows
    for j = 1 : rows
        z(j,:) = z(j,:) / (sigma(j) + a_tol);
    end
elseif length(sigma) == cols
    for j = 1 : cols
        z(:,j) = z(:,j) / (sigma(j) + a_tol);
    end
end

end

