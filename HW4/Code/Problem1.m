% Shabab Ahmed 
% ECON 525 HW 4
%--------------------------------------------------------------------------
% Problem 1(b): two step
%--------------------------------------------------------------------------
%%

% Loading probability of i_t given x_t
load Prob.mat

% Set seed 
rng(100)

% Simulating new data 
n_sims = 500000;
x_data = zeros(n_sims, 1); 
i_data = zeros(n_sims, 1);
x_data(1) = 0; 

% Lambda 
lam = 0.8;

for i=2:n_sims 
    val = binornd(1, lam); 
    if i_data(i-1) == 1
        x_data(i) = 0+val;
    else 
        x_data(i) = min(x_data(i-1)+val, 10); 
    end
    i_data(i)= binornd(1, P_ivf(x_data(i)+1));
end 

simulated_data = [x_data, i_data];

T= length(x_data); 

% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);

% Probability of replacement at each state 
P_hat = zeros(11,1); 

for i=1:nx
    I = find(x_data == x_grid(i)); % sum of 1(x_t=x)
    denom = length(I); % denominator in p_hat
    num = length(find(i_data(I)==1)); % numerator in p_hat
    P_hat(i) = num/denom; 
    % To ensure we can take logs
    if P_hat(i) == 1
        P_hat(i) = 0.99;
    end 
end



% Calculation of p_hat

x_p = x_data;
i_p = i_data; 
x_p(T)=[];
i_p(T)=[];

p_hat_0 = zeros(11,11); 
p_hat_1 = zeros(11,11);

for i=1:nx
    for j=1:nx
        I = find(x_p == x_grid(i)); 
        i_1 = find(i_p(I)==1);
        i_0 = find(i_p(I)==0); 
        denom_0 = length(i_0); 
        denom_1 = length(i_1);
        xi_1 = I(i_1);
        xi_0 = I(i_0); 
        xi_1(xi_1 == 499999) = [];
        xi_0(xi_0 == 499999) = [];
        x_prime_1 = find(x_p(xi_1+1)==x_grid(j));
        x_prime_0 = find(x_p(xi_0+1)==x_grid(j));
        num_1 = length(x_prime_1); 
        num_0 = length(x_prime_0);
        p_hat_0(i,j) = num_0/denom_0; 
        p_hat_1(i,j) = num_1/denom_1;
    end 
end 

%%
%--------------------------------------------------------------------------
% Problem 1(b)
%--------------------------------------------------------------------------
first = 0;
last =20; 
int = 1; 

% Values of parameters to be searched over 

theta1 = (first:int:last).*(1/10);
theta2= (first:int:last).*(1/100);
theta3 = (first:int:last).*(5/10);
len_theta = length(theta1);
% Creating list of parameters to be searched over 

params = [];
for i=1:len_theta
    for j = 1:len_theta 
        for k=1:len_theta
            theta = [theta1(i), theta2(j), theta3(k)];
            params = [params; theta];
        end
    end
end 

% Euler constant
u = 0.5772;

% Creating G 
G = zeros(nx, nx);
for i=1:nx
    P_x = P_hat(i); 
    for j=1:nx
        G(i,j)= P_x*(p_hat_1(i,j))+(1-P_x)*(p_hat_0(i,j));
    end
end


len_params = length(params);

% saving likelihood for different paramenters
l_params = []; 

% Grid search 
for i = 1:len_params
    parameter = params(i,:);
    theta_1 = parameter(1);
    theta_2 = parameter(2);
    theta_3 = parameter(3); 
    [V_0, V_1] = IVF(P_hat, p_hat_0, p_hat_1,G, theta_1, theta_2, theta_3);
    % calculating likelihood
    L_theta = Likelihood(x_data, i_data, p_hat_0, p_hat_1, V_0, V_1); 
    l_params = [l_params; L_theta];
end

[L, ind] = max(l_params);

% Parameters that maximize likelihood

param_MLE = params(ind,:);

%%
%--------------------------------------------------------------------------
% Problem 1(c): forward simulation
%--------------------------------------------------------------------------
rng(100)
diff_V_barhat = log(P_hat)-log(1-P_hat);
S = 100; 
T=200; 

% Simulating the data 
data = []; 
data_x = []; 
data_i = [];
eps_0 = []; 
eps_1 = [];
for i=1:11
    disp(i)
    x= x_grid(i);
    I_data=[];
    X_data = [];
    for j = 1:S
        I_t = zeros(T,1);
        x_t = x*ones(T,1);
        for t=1:(T)
            lam = random('Uniform',0,1);
            epsilon_0 = -random('ev',0,1); 
            epsilon_1 = -random('ev',0,1); 
            eps_0 = [eps_0; epsilon_0]; 
            eps_1 = [eps_1;epsilon_1];
            diff = epsilon_0-epsilon_1; 
            if diff_V_barhat(x_t(t)+1)>diff
                I_t(t) = 1; 
            end 
            if t < T
                if I_t(t) == 1 
                    if lam > p_hat_1(x_t(t)+1,1)
                        x_t(t+1) = 1; 
                    else 
                        x_t(t+1) = 0; 
                    end
                else 
                    if lam > p_hat_0(x_t(t)+1,x_t(t)+1)
                        x_t(t+1) = x_t(t)+1; 
                    else 
                        x_t(t+1) = x_t(t); 
                    end 
                end
            end
        end 
        X_data = [X_data;x_t];
        I_data = [I_data; I_t];    
      
    end 
    data_x = [data_x; X_data];
    data_i = [data_i; I_data];
end
data = [data_x, data_i];


%%
% Calculating value functions 

A_vec= zeros(11,1); 
B_vec = zeros(11,1);
C_vec = zeros(11,1);
D_vec = zeros(11,1);
beta= 0.95;

for i=1:length(x_grid) 
    len1 = length(data_x)/11*(i-1)+1;
    len2 = length(data_x)/11*i;
    A = 0; 
    B = 0; 
    C = 0;
    D = 0;
    for j = len1:len2
        t = mod(j, T); 
        if t==0
            t = T;
        end
        b = beta^t; 
        B = B+b*data_i(j); 
        C = C+b*(eps_1(j)*data_i(j)+eps_0(j)*(1-data_i(j)));
        A = A +b*data_x(j)*(1-data_i(j));
        D = D + b*(data_x(j)*(1-data_i(j)))^2; 
    end 
    A_vec(i) = (1/S)*A; 
    B_vec(i) = (1/S)*B;
    C_vec(i) = (1/S)*C;
    D_vec(i) = (1/S)*D; 
end

%%
% MLE 

first = 0;
last =20; 
int = 1; 

% Values of parameters to be searched over 

theta1 = (first:int:last).*(1/10);
theta2= (first:int:last).*(1/100);
theta3 = (first:int:last).*(4/10);
len_theta = length(theta1);
% Creating list of parameters to be searched over 

params = [];
for i=1:len_theta
    for j = 1:len_theta 
        for k=1:len_theta
            theta = [theta1(i), theta2(j), theta3(k)];
            params = [params; theta];
        end
    end
end 

len_params = length(params);

l_params_fs = []; 

% Grid search 
for i = 1:len_params
    disp(i)
    parameter = params(i,:);
    theta_1 = parameter(1);
    theta_2 = parameter(2);
    theta_3 = parameter(3); 
    cont_payoff = -theta_1*A_vec-theta_2*D_vec-theta_3*B_vec+C_vec;
    V_0 = -theta_1*x_grid-theta_2*x_grid.^2+p_hat_0*cont_payoff; 
    V_1 = -theta_3+p_hat_1*cont_payoff;
    % calculating likelihood
    L_theta_fs = Likelihood(x_data, i_data, p_hat_0, p_hat_1, V_0, V_1); 
    l_params_fs = [l_params_fs; L_theta_fs];
end

[L_fs, ind_fs] = max(l_params_fs);

% Parameters that maximize likelihood
param_MLE_fs = params(ind_fs,:);
%%
% Calculating value function solutions 
theta_1 = param_MLE_fs(1); 
theta_2 = param_MLE_fs(2);
theta_3 = param_MLE_fs(3);
cont_payoff = -theta_1*A_vec-theta_3*B_vec+C_vec;
V_0 = -theta_1*x_grid-theta_2*x_grid.^2+p_hat_0*cont_payoff
V_1 = -theta_3+p_hat_1*cont_payoff