function [r2] = rSquared(Y, F, scal_crit)
%% Description
% Computes the R-squared values.
% INPUT
% - Y: (samples x variables) matrix of original data.
% - F: (samples x variables) matrix of predicted data.
% - scal_crit: scaling criterion (0 = no scaling, only centering)

%% Input
a_tol = 1e-16;
% Check they have the same size
[m, n] = size(Y);
if m ~= size(F,1) || n ~= size(F,2)
    error('Input matrices must have the same size.');
end
% Center-scale data
if nargin < 3
    scal_crit = 0; % Subtract the mean only
end
Y = center_scale_data(Y, scal_crit);
F = center_scale_data(F, scal_crit);

%% Main
SStot = sum(Y.^2, 1);
SSres = sum((Y - F).^2, 1);

%% Output
r2 = 1 - SSres./(SStot + a_tol);

end

