function [y] = removeHoles(x, t)
%% Description
% This function removes the missing number in a series of integers.
% Example: [1 2 4 5] becomes [1 2 3 4].
% Input:
% - x: (N x 1) array of integers
% - y: (N x 1) array of integers
%

%% Main
u = unique(x);
n_ints = length(u); % Number of unique elements in x
% Loop over values in u
for ii = n_ints - 1 : -1 : 1
    % u(ii) must equal u(ii+1) - 1
    if u(ii) + 1 ~= u(ii+1)
        m = u(ii+1) - u(ii); % Hole's size
        I = (x > u(ii)); % Values in x who are greater than u(ii)
        x(I) = (x(I) - m) + 1; % Get rid of the hole (but add 1)
    end
end
% min(y) must be 1
if nargin < 2
    t = false;
end
if t
    m = (min(x) - 1);
    x = x - m;
end

%% Output
y = x;

end

