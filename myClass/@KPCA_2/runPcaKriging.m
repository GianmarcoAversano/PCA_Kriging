function runPcaKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end
% Kriging
n_args = length(varargin);
if n_args > 0 && ~isempty(varargin{1})
    % User-supplied kriged scores
    obj.kriged_pca_scores = varargin{1};
else
    % Kriging on the PCA scores
    fprintf('Kriging on the PCA scores is in progress...\n');
    temp = 'pca';
    test_data = [];
    if ~isempty(obj.original_data)
        test_data = obj.pca_modes' * obj.original_data;
    end
    % Empty properties
    obj.pcaKrigingModel = {};
    obj.mse = [];
    % Train model
    evalc('obj.update_Kriging(obj.pca_scores, temp, test_data);');
    fprintf('Kriging on the PCA scores has terminated.\n');
end
% Get the predictions
obj.getKpcaPredictions();
% Terminate here if this is a LPCA object
if obj.is_local
    return
end
% Estimate the errors
obj.getKpcaPredictionsErrors();

end


%% Obsolete pieces

% TARGET BY TARGET 
%     if obj.targetBYtarget
%         k = cell(size(obj.pca_scores,1), 1);
%         mse = zeros(size(obj.pca_scores,1), size(obj.prediction_points,1));
%         for i = 1 : size(obj.pca_scores,1)
%             if ~isempty(obj.prediction_points);
%                 [obj.kriged_pca_scores(i,:), k{i}, mse(i,:)] = obj.update_Kriging(obj.pca_scores(i,:), temp, test_data(i,:));
%             else
%                 obj.update_Kriging(obj.pca_scores(i,:), temp, test_data(i,:));
%             end
%         end
%         obj.pcaKrigingModel = k;
%         obj.pcaKrigingMSE = mse;
%     else
%         obj.update_Kriging(obj.pca_scores, temp, test_data);
%     end

