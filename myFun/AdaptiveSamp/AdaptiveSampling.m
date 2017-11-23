function [varargout] = AdaptiveSampling(starting_points, starting_values, cdt_points, varargin)
%% Description
% Adaptive Sampling for basis improvements.
% INPUTS
%     starting_values: var x obs
%     cdt_points: obs x var
%     starting_points: obs x var
% OUTPUTS
%     point: new point to add
%

%% Input check
% varargin size
n_args = length(varargin); in_count = 1;
% varargin check
reg_based = false;
if n_args > in_count
    in_count = in_count + 1;
    reg_based = varargin{in_count};
end
% Be sure that no value in p0 is in p
I = ismember(cdt_points, starting_points, 'rows');
cdt_points(I) = [];
% Number of points
N0 = size(starting_points,1);
N_cdt = size(cdt_points,1);

%% Relative Influence
if reg_based
    % Kriging-based
    trendFun = varargin{in_count}; in_count = in_count + 1;
    corrFun = varargin{in_count}; in_count = in_count + 1;
    is_slow = false;
    Infl_rel = samplingInfluence_Kriging(starting_points, starting_values', trendFun, corrFun, [], is_slow);
else
    % PCA-based
    Infl_rel = samplingInfluence_PCA(starting_points, starting_values);
end

%% Enrichment Potential
% Memory allocation
Pot = zeros(N_cdt, 1);
dist = zeros(N0, 1);
% Potential of each candidate point
for i = 1 : N_cdt
    % Evaluate the distance of all starting points from this candidate
    % points
    this_point = cdt_points(i,:);
    parfor j = 1 : N0
        dist(j) = norm(this_point - starting_points(j,:));
    end
    % Find the index of the closest point (the starting point that is the
    % closest to this candidate point
    [~, minloc] = min(dist);
    % Choose the j for which it is max
    Pot(i) = dist(minloc) * Infl_rel(minloc);    
end
% Candidate point to be chosen: get the index of the candidate point that 
% maximizes the Pot
[~, maxloc] = max(Pot);

%% Output
out_count = 0;
if nargout > out_count
    out_count = out_count + 1;
    varargout{out_count} = maxloc;
end
if nargout > out_count
    out_count = out_count + 1;
    varargout{out_count} = Pot;
end
if nargout > out_count
    out_count = out_count + 1;
    varargout{out_count} = Infl_rel;
end

end



