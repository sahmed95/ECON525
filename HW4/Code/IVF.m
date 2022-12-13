function [V_0, V_1] = IVF(P_hat, p_hat_0, p_hat_1, G,theta_1, theta_2, theta_3)
% State space
x_min = 0;
x_max = 10;
interval = 1;
x_grid = (x_min:interval:x_max)';
nx = length(x_grid);


% Utility functions
U_0 = @(x)-theta_1*x-theta_2*x.^2;
U_1 = -theta_3;

% Parameter values 
beta = 0.95; 

% Euler constant
u = 0.5772;

eps0_hat = u-log(1-P_hat);
eps1_hat = u-log(P_hat);

% Creating Pi matrix 

Pi = zeros(nx,1);
for i=1:nx
    P_x = P_hat(i); 
    eps0_x = eps0_hat(i);
    eps1_x = eps1_hat(i);
    x = x_grid(i);
    Pi(i) = P_x*(U_1+eps1_x)+(1-P_x)*(U_0(x)+eps0_x);
end


% calculating new value function
I = eye(nx,nx);
V = (I-beta*G)\Pi;
V_0 = U_0(x_grid)+beta*p_hat_0*V; 
V_1 = U_1+beta*p_hat_1*V;



