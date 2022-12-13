% Shabab Ahmed 
% ECON 525 HW 2-3
%--------------------------------------------------------------------------
% Problem 2 (c)
%--------------------------------------------------------------------------
%%
% Loading draw
load draw.csv.csv 
lambda = draw_csv(:,2);
i_draw = draw_csv(:,1);

% Loading probability 
load Prob.mat

% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);

% Time period
T = 5000; 

% Storing simulated data with intial values of zero
x_data = zeros(T,1);
i_data = zeros(T,1);


bern = zeros(T, 1);
for i=1:T
    if lambda(i)>0.2
        bern(i) = 1;
    end 
end 
if i_draw(1) > 1-P_ivf(x_data(1)+1)
    i_data(1) = 1;
end 

% Simulating the data
for i=2:T
    if i_data(i-1) == 1
        if bern(i) == 1
            x_data(i) = 1; 
        else 
            x_data(i) = 0;
        end 
    else
        if bern(i) ==1
            x_data(i) = min([x_data(i-1)+1,10]); 
            
        else 
            x_data(i) = x_data(i-1);
          
        end
    end 
    if i_draw(i) > 1-P_ivf(x_data(i)+1)
            i_data(i) = 1;
    end
end 

simulated_data = [x_data, i_data];

% Saving data for HW 4

save('data.mat', "simulated_data")

% Removing observations with x_t = 10 and i_t = 0
x_10 = x_data; 
i_10 = i_data;
ind_10 = find(x_data == 10);
if i_data(ind_10) == 0
    x_10(ind_10) =[]; 
    i_10(ind_10) = [];
end 

sim_10 = [x_10, i_10];

% Calculating lambda using MLE
len = length(sim_10);
lam_sum = 0;
for i=2:len
    val = i_10(i-1)*x_10(i)+(1-i_10(i-1))*(x_10(i)-x_10(i-1));
    lam_sum = lam_sum +val;
end

lambda = (1/len)*lam_sum;

save('lambda.mat', 'lambda')
%%
%--------------------------------------------------------------------------
% Problem 2 (d)
%--------------------------------------------------------------------------
first = 0;
last =10; 
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

len_params = length(params);

% saving likelihood for different paramenters

l_params = []; 

% Grid search 

for i = 1:len_params
    parameter = params(i,:);
    theta_1 = parameter(1);
    theta_2 = parameter(2);
    theta_3 = parameter(3); 
    [Pi_0, Pi_1, P] = IVF(lambda, theta_1, theta_2, theta_3); 
    % calculating likelihood
    L_theta = Likelihood(x_data, i_data, Pi_0, Pi_1, P); 
    l_params = [l_params; L_theta];
end

[L, ind] = max(l_params);

% Parameters that maximize likelihood

param_MLE = params(ind,:);

save('estim_param.mat', "param_MLE")


%%


