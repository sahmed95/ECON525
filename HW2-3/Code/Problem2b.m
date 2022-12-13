% Shabab Ahmed 
% ECON 525 HW 2-3

%%
%--------------------------------------------------------------------------
% Problem 2 (b) (1): Choice-specific value function formulation
%--------------------------------------------------------------------------
% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);

% Parameter values 
theta_1 = 0.3;
theta_2 = 0.0;
theta_3 = 4.0;
lambda = 0.8; 
beta = 0.95; 


% Transition matrix for i_t = 0

Pi_0 = zeros(11,11); 
for i = 1:10
    Pi_0(i,i) = 1-lambda; 
end 

for i=1:10
    Pi_0(i, i+1) = lambda;
end

% if x_t = 10, we stay there 

Pi_0(11,11)= 1;

% Transition matrix for i_t = 1 

Pi_1 = zeros(11,11);

for i = 1:nx
    Pi_1(i, 1) = 1-lambda; 
    Pi_1(i, 2) = lambda; 
end
% Convergence criterion 
crit0 = 1; 
crit1 = 1;
tol = 1e-5;
max_iter = 10000; 

% Euler constant 
u = 0.5772;

v_0 = zeros(nx,1); % value function for i_t = 0; 
v_1 = zeros(nx,1); % value function for i_t = 1; 
tv_0 = zeros(nx,1);
tv_1 = zeros(nx,1);
tv = zeros(nx,1);
v = zeros(nx,1);
dr = zeros(nx,1); % decision rule (will contain indices)
count=0; % to keep track of iterations
while crit0>tol && crit1>tol && count < max_iter  
    count=count+1;
    disp([num2str(count), 'th iteration ']);
    val = u+log(exp(v_0)+exp(v_1)); 
    for i=1:nx
        % Utility functions
        U_0 = -theta_1*x_grid(i)-theta_2*x_grid(i)^2;
        U_1 = -theta_3;
        
        % Value functions 
        tv_0(i)= U_0 + beta*Pi_0(i,:)*val; 
        tv_1(i) = U_1+ beta*Pi_1(i,:)*val;
        tv(i) = max([tv_0(i); tv_1(i)]);
    end  
    crit0 = max(abs(tv_0-v_0)); 
    crit1 = max(abs(tv_1-v_1)); 
    v_0 = tv_0; % Update the value function
    v_1 = tv_1; 
    v = tv; 
end 


% Replacement probabilities 
P = zeros(nx,1); 
for i=1:nx
    num = exp(v_1(i));
    denom = exp(v_0(i))+num;
    P(i) = num/denom; 
end 


%%

%--------------------------------------------------------------------------
% Problem 2 (b) (2): Integrated value function formulation
%--------------------------------------------------------------------------
% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);

% Parameter values 
theta_1 = 0.3;
theta_2 = 0.0;
theta_3 = 4.0;
lambda = 0.8; 
beta = 0.95; 

% Euler constant
u = 0.5772;

% Transition matrix for i_t = 0

Pi_0 = zeros(11,11); 
for i = 1:10
    Pi_0(i,i) = 1-lambda; 
end 

for i=1:10
    Pi_0(i, i+1) = lambda;
end

% if x_t = 10, we stay there 

Pi_0(11,11)= 1;

% Transition matrix for i_t = 1 

Pi_1 = zeros(11,11);

for i = 1:nx
    Pi_1(i, 1) = 1-lambda; 
    Pi_1(i, 2) = lambda; 
end

% Convergence criterion 
crit = 1; 
tol = 1e-5;
max_iter = 10000; 


% Utility functions
U_0 = @(x)-theta_1*x-theta_2*x.^2;
U_1 = -theta_3;

% Initial replacement probability and value function
P_ivf = ones(nx,1)*0.5; 
v_ivf = zeros(nx,1); 
v_0_ivf = zeros(nx,1);
v_1_ivf = zeros(nx,1);
tv_ivf = zeros(nx,1);
count = 0;
while crit>tol && count < max_iter  
    count = count+1;
    disp([num2str(count), 'th iteration ']);
    
    eps0_hat = u-log(1-P_ivf);
    eps1_hat = u-log(P_ivf);

% Pi and G
    Pi = zeros(nx,1);
    G = zeros(nx, nx);
    for i=1:nx
        P_x = P_ivf(i); 
        eps0_x = eps0_hat(i);
        eps1_x = eps1_hat(i);
        x = x_grid(i);
        Pi(i) = P_x*(U_1+eps1_x)+(1-P_x)*(U_0(x)+eps0_x);
        for j=1:nx
            G(i,j)= P_x*(Pi_1(1,j))+(1-P_x)*(Pi_0(i,j));
        end
    end
    % calculating new value function
    I = eye(nx,nx);
    tv_ivf = (I-beta*G)\Pi;
    crit = max(abs(tv_ivf-v_ivf));
    % updating V_0 and V_1
    v_0_ivf = U_0(x_grid)+beta*Pi_0*tv_ivf; 
    v_1_ivf = U_1+beta*Pi_1*tv_ivf;
    v_ivf = tv_ivf;
    % updating P
    P_ivf = (exp(v_1_ivf))./(exp(v_1_ivf)+exp(v_0_ivf));
end 

save("Prob.mat", 'P_ivf')
%%

