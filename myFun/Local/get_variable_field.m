function [Y_var] = get_variable_field(Y, mesh, vars, variable, X, x, varargin)
%% Description:
% Y: obs x var 

%% Main
if ischar(variable)
    [Y_var] = get_variable_field_local(Y, mesh, vars, variable, X, x);
elseif iscell(variable)
    Y_var = [];
    l = length(variable);
    for ii = 1 : l
        [temp] = get_variable_field_local(Y, mesh, vars, variable{ii}, X, x);
        Y_var = [Y_var, temp];
    end
else
    error('There seems to be a problem with the input VARIABLE.');
end
end

function [Y_var] = get_variable_field_local(Y, mesh, vars, variable, X, x, varargin)
% Y: obs x var 
n_mesh = size(mesh,1);
n_vars = length(vars);
rows = size(Y,1);
for ii = 1 : n_vars
    if strcmp(variable, vars{ii})
        idx = ii;
        break
    end
end
x1 = (idx - 1) * n_mesh + 1;
x2 = idx * n_mesh;
Y_var = Y(:,x1:x2);
if nargin < 5
    return
end
I = ismember(X, x, 'rows');
Y_var = Y_var(I,:);
end



