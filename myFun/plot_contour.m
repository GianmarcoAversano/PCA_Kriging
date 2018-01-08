function [varargout] = plot_contour(mesh, y_data, preproc_flag, new_fig, flipped, varargin)
% DESCRIPTION: Works only for Zhy Li's case.

% Input
cont = false;
if nargin < 3 || isempty(preproc_flag)
    preproc_flag = false;
end
if nargin < 4 || isempty(new_fig)
    new_fig = true;
else
    new_fig = false;
end
if nargin < 5
    flipped = false;
end
% Make mesh denser
if preproc_flag
    % Get unique coordinates
    x = (unique(mesh(:,1)));
    y = (unique(mesh(:,2)));
    % Make the mesh more dense 
    deltax = abs(x(2) - x(1));
    deltay = abs(y(2) - y(1));
    delta = min([deltax, deltay]) / 3.0;
    x = mesh(1,1) : delta : mesh(end,1);
    y = mesh(1,2) : delta : mesh(end,2);
else
    % Get unique coordinates
    x = (sort(mesh(:,1)));
    y = (sort(mesh(:,2)));
end
% Plot original data
[xx, yy] = meshgrid(x, y);
F = scatteredInterpolant(mesh(:,1), mesh(:,2), y_data(:));
zz = F(xx, yy);
if flipped
    F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
    zz2 = F(xx, yy);
    zz = [flip(zz,2), zz2];
    [xx, yy] = meshgrid([flip(-x); x], y);
end
if new_fig
    figure();
end
if cont
    contourf(xx, yy, zz);
else
    if flipped
        imshow(zz, [], 'XData', x, 'YData', y, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on;
    else
        imshow(zz, [], 'XData', [flip(-x); x], 'YData', y, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on; 
    end
end
xlabel('x [m]');
ylabel('y [m]');
end



% Make the mesh more dense 
% deltax = abs(x(2) - x(1));
% deltay = abs(y(2) - y(1));
% delta = min([deltax, deltay]) / 3.0;
% x = min((obj.mesh(:,1))) : delta : max((obj.mesh(:,1))); x = x(:);
% y = min((obj.mesh(:,2))) : delta : max((obj.mesh(:,2))); y = y(:);
% % Get data to plot
% [xx, yy] = meshgrid(x, y);
% F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y2plot, method);
% zz = F(xx, yy);
% F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
% zz2 = F(xx, yy);
% zz = [flip(zz,2), zz2];
% [xx, yy] = meshgrid([flip(-x); x], y);
% % Plot
% figure();
% if cont
%     y_min = min(y2plot); y_max = max(y2plot); % Range
%     v = y_min + (0 : .05 : 1) * (y_max - y_min); % Color values
%     [~, h, ~] = contourf(xx, yy, zz, v); 
%     set(h(:), 'LineStyle', 'none');
% else
%     imshow(zz, [], 'XData', [flip(-x); x], 'YData', y, 'Colormap', jet);
%     axis xy; % Flip the figure
%     axis on; 
% end
