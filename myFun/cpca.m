function [gamma, debug] = cpca(X0, my_constraint, scores, modes, x_ave, x_sig, cpca_guess_corrector)
%% Description
% X0: (n_obs x n_var) data from which MODES and SCORES were extracted, which PCA was
%                    performed on.
% my_constraint: function handle of the non linear constraint to be
%                respected.
% scores: (n_obs x n_pc) PCA scores
% modes: PCA modes
% x_ave: vector of mean values
% x_sig: vector of scaling factors
%

%% Input
if isempty(my_constraint)
    my_constraint = @cpca_con;
end
if nargin < 7
    cpca_guess_corrector = 1;
end

%% Evaluating the CPCA scores:
% Starting points (or guesses)
gamma0 = cpca_guess_corrector * scores; 
disp('--- Evaluating CPCA scores ---');
% Solver:
if isa(my_constraint, 'function_handle') 
    % Non-linear constrained optimization problem
    gamma = solv(my_constraint, scores', modes, [], [], gamma0', X0', x_ave, x_sig);
elseif isa(my_constraint, 'cell') || isa(my_constraint, 'double')
    % Linear constrained problem
    [gamma, debug] = solv2(obj, gamma0);
end

%% Output:
if nargout > 1
    varargout{1} = debug;
end
end


%------- Solving for the gamma's ---------------%
function [gamma] = solv(constraint_fun, scores, modes, lb, ub, gamma0, X, x_ave, x_sig)
n_obs = size(X, 2);
k = size(scores, 1);
% Evaluation of each gamma(:,i)
gamma = zeros(k, n_obs); % Memory allocation
% Possible algorithms for the solver of the constrained non-linear problem:
x = {'interior-point', 'sqp'};  % {'interior-point'; 'sqp'; 'active-set'};
% Loop repetitions
parfor i = 1 : n_obs       
    % Starting guess:
    guess = gamma0(:,i);
    % Solving
    j = 1;
    goforth = true;
    while goforth && (j <= length(x))
        % Options:
        options = setOptions(x{j});
        % Use "FMINCON":
        [gammasol, fval, exitflag] = fmincon(...
            @(gamma) ObjFun(gamma, X(:,i), modes),... % Objective function
            guess,...       % Starting guess
            [],[],...       % Linear inequality costraints
            [],[],...       % Linear equality costraints 
            lb,ub,...       % Lower and upper bounds (lb, ub)
            @(gamma) constraint_fun(gamma, modes, x_ave, x_sig),... % Non-linear function for costraints
            options);       % Options
        guess = gammasol;
        % Exitflag
        if exitflag < 1
            j = j + 1;
        else
            goforth = false;
        end 
    end
    % Solution:
    gamma(:,i) = gammasol; 
    % Show real time to prove work is in progress
    fprintf('Solved for case %i out of %i\n', i, n_obs);
end % parfor
end % end of solv()
%------- End: Solving for the gamma's ---------------% 


%------- Set options ---------------%
function options = setOptions(s)
switch s
    case 'interior-point'
        MaxIter = 5e5;     TolFun = 1e-12;     TolX = 1e-12; 
        TolCon = 1e-8; 
    case 'sqp'
        MaxIter = 1e6;     TolFun = 1e-15;     TolX = 1e-15;
        TolCon = 1e-12;
    case 'active-set'
        MaxIter = 5e2;     TolFun = 1e-18;     TolX = 1e-18;
        TolCon = 1e-12;
end
MaxFunEvals = 1e11;
options = optimoptions('fmincon',...
        'Algorithm',s,...
        'MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter, 'TolFun',TolFun,...
        'TolX',TolX, 'TolCon',TolCon,...
        'Display','off'); % 'iter', 'final-detailed'
end
%------- End: Set options ---------------%


%------- ObjFun ----------------------%
function f = ObjFun(gamma, y, modes) 
y_rec = modes * gamma;
err = y - y_rec;
h = norm(y) + eps;
f = .5 * norm( err ) / ( h );
end
%------- End: ObjFun ---------------%

%------- cpca_con ------------------%
function c = cpca_con(cpca_scores, modes, x_ave, x_sig)
%% Description:
% This is an non-linear constraints which requires the rebuilt data vector
% y to have all positive values. It can be a common constraint to be used,
% so it is present here by default.
%
% OUTPUTS:
% c(x) is the array of nonlinear inequality constraints at x: 
% fmincon attempts to satisfy c(x) <= 0 for all entries of c.
%
% ceq(x) is the array of nonlinear equality constraints at x. 
% fmincon attempts to satisfy ceq(x) = 0 for all entries of ceq.
% 
% For example, x = fmincon(@myfun, x0, A, b, Aeq, beq, lb, ub, @mycon)
%
% INPUTS:
% cpca_scores: vector of CPCA scores
% modes: PCA modes
% x_ave: mean values (vector)
% x_sig: scaling factors (vector)
%

%% Pre-proc:
a_tol = 1e-16;
x_ave(x_ave < a_tol) = 0; % To avoid problems
% Column vectors
cpca_scores = cpca_scores(:);
x_ave = x_ave(:);
x_sig = x_sig(:);

%% Non-linear inequalities:
% The vector y is rebuilt
y = modes * cpca_scores;
y = x_ave + x_sig .* y; 
% All values have to be > 0
c = - min(y);                       

%% Non-linear equalities:
ceq = []; % No non-linear inequality constraints.

end
%------- cpca_con ------------------%




