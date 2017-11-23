function [varargout] = getError(Y_ref, Y_app, dim, sc, varargin)
%% Description
% Input matrices Y_ref and Y_app are to be interpreted as (variables \times
% observations).
%

%% Input
if nargin < 3
    dim = 0;
end
if nargin < 4
    sc = 0;
end

%% Main
a_tol = 1e-9;
% Get scaling factors
[~, ~, sf] = center_scale_data(Y_ref, sc); % 1: std; 2: range
% Evaluate difference matrix
ER2 = abs(Y_ref - Y_app);
% Scale the difference matrix
ER2 = ER2 .* repmat(sf, 1, size(Y_ref,2)); % Multiply
ER2 = ER2 ./ (repmat(sf.^2, 1, size(Y_ref,2)) + a_tol); % Divide
% Get a vector
ER1 = mean(ER2, 2);
% Get a scalar
ER0 = mean(ER1);

%% Output
if dim == 0
    varargout{1} = ER0;
elseif dim == 1
    varargout{1} = ER1;
elseif dim == 2
    varargout{1} = ER2;
end

end

