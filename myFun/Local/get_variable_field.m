function [Y_var] = get_variable_field(Y, mesh, vars, variable, X, x, varargin)
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




