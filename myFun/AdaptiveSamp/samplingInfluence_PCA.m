function [Infl_rel, varargout] = samplingInfluence_PCA(X, Y, cent_crit, scal_crit, varargin)
%% Description
%{
part of: "Adaptive Sampling for basis improvements"

Evaluation of the Influence of each snapshot: 
- on the PCA basis;
- on the regression model;

INPUTS
    Y: matrix of observations (column-vectors)
    points: input space locations of the observations

OUTPUTS
    Infl_rel: relative influence of each observation
%}

%% Input check
% Centergin criterion
if nargin < 3
    cent_crit = 1;
end
% Scaling criterion
if nargin < 4
    scal_crit = 1;
end
% Number of points
n_points = size(X,1);

%% PCA 
pca_model = pcaFun(Y, cent_crit, scal_crit); % Create a PCA structure

%% Snapshot influences on the modes
a_tol = 1e-16;
% Memory allocation
pod_nosnap = cell(n_points,1);
Infl = [];
% Create the matrix of snapshot influences on the modes
for j = 1 : n_points
    fprintf('[samplingInfluence_PCA] %i out of %i \n', j, n_points);
    % Create the matrix with the missing snapshot
    Yj = Y;   
    Yj(:,j) = []; 
    % PCA with a missing snapshot
    pod_nosnap{j} = pcaFun(Yj, cent_crit, scal_crit); 
    % Infl(i,j): influence of snapshot j on the pca mode i
    for i = 1 : pod_nosnap{j}.k
        modeprod = abs( pca_model.modes(:,i)' * pod_nosnap{j}.modes(:,i) );
        Infl(i,j) = 1 / (modeprod + a_tol) - 1; 
    end
end

%% Influence of the snapshots on the modal basis
% Memory allocation
Infl_s = zeros(length(pod_nosnap), 1);
% Influence of the snapshot j on the modal basis
for j = 1 : length(pod_nosnap)
    sv = (pod_nosnap{j}.eigenv).^(.5); sv = sv(:); % Singular value
    Infl_s(j,1) = sv' * Infl(:,j);
end

%% Output
% Relative influence of the jth snapshot on the modal basis
Infl_rel = Infl_s / norm(Infl_s, 1);
% Return output
if nargout > 1 
    varargout{1} = Infl;
end


end

function this = pcaFun(Y, cent_crit, scal_crit)
Y = Y';
[this.Ycs, this.m] = center(Y, cent_crit); % Y_ave is a vector
[this.Ycs, this.d] = scale(this.Ycs, Y, scal_crit); % Y_gamma is a vector
this.Ycs = this.Ycs';
% Apply PCA (Using Matlab's default toolbox)
[coeff, scores, latent] = pca(this.Ycs');
% Approximation order
this.k = size(coeff, 2);
% Modes
this.modes = coeff;        
% Scores
this.a = scores';           
% Eigenvalues
this.eigenv = latent;       
end



