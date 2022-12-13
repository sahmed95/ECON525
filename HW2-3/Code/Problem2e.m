% Shabab Ahmed 
% ECON 525 HW 2-3
%--------------------------------------------------------------------------
% Problem 2 (e) (i)
%--------------------------------------------------------------------------
%%
% Setting seed
rng(1000); 

% Loading estimated parameters
load estim_param.mat
theta_1 = param_MLE(1);
theta_2 = param_MLE(2); 
theta_3 = param_MLE(3);

% Lambda from MLE
load lambda.mat


% Number of data points
T = 5000;

% Simulating new data 
x_new = zeros(T,1);
i_new = zeros(T,1);

% Solving DP problem 

[Pi_0, Pi_1, P] = IVF(lambda,theta_1, theta_2, theta_3);

% First i_t value 
i_new(1) = binornd(1, P(1)); 

for i=2:T
    if i_new(i-1) == 1
        x_new(i) = binornd(1, lambda); 
    else
        lambda_new = binornd(1,lambda);
        x_new(i) = min(x_new(i-1)+lambda_new, 10); 
    end 
    i_new(i)= binornd(1, P(x_new(i)+1)); 
end 

prob_replacement_sim = (1/T)*sum(i_new);
%%
%--------------------------------------------------------------------------
% Problem 2 (e) (ii): using steady state
%--------------------------------------------------------------------------
% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);

% Euler constant 
u = 0.5772;

eps_1 = u-log(P); 
eps_0 = u-log(1-P);

% Utility functions 
U_0 = @(x)-theta_1*x-theta_2*x.^2;
U_1 = -theta_3;


Pi = zeros(nx,1);
G = zeros(nx, nx);
for i=1:nx
    P_x = P(i); 
    eps0_x = eps_0(i);
    eps1_x = eps_1(i);
    x = x_grid(i);
    Pi(i) = P_x*(U_1+eps1_x)+(1-P_x)*(U_0(x)+eps0_x);
    for j=1:nx
        G(i,j)= P_x*(Pi_1(1,j))+(1-P_x)*(Pi_0(i,j));
    end
end
% Calculating steady state

[V, D] = eig(G');                          
ii = find(abs(diag(D) - 1) < 1E-8, 1);      % Find first unit eigenvalue           
A = V(:,ii) / sum(V(:,ii));                 % Normalize unit eigenvector
assert(max(abs(A' - A'*G)) < 1E-12)        % Verify dist. is stationary
if sum(abs(diag(D) - 1) < 1E-8) > 1         % Check for uniqueness
   warning('A not unique')
end

prob_replacement_steady_state = dot(A,P);
%%
%--------------------------------------------------------------------------
% Problem 2 (e) (iii): using simulation
%--------------------------------------------------------------------------
rng(1000)
% New theta_3
theta_3 = 0.9*param_MLE(3);

% Number of data points
T = 5000;

% Simulating new data 
x_new_subsidy = zeros(T,1);
i_new_subsidy= zeros(T,1);

% Solving DP problem 

[Pi_0, Pi_1, P_subsidy] = IVF(lambda,theta_1, theta_2, theta_3);

% First i_t value 
i_new_subsidy(1) = binornd(1, P(1)); 

for i=2:T
    if i_new(i-1) == 1
        x_new_subsidy(i) = binornd(1, lambda); 
    else
        lambda_new = binornd(1,lambda);
        x_new_subsidy(i) = min(x_new_subsidy(i-1)+lambda_new, 10); 
    end 
    i_new_subsidy(i)= binornd(1, P_subsidy(x_new_subsidy(i)+1)); 
end 

prob_replacement_sim_subsidy = (1/T)*sum(i_new_subsidy);

%%

