% Shabab Ahmed 
% ECON 525 HW 5
%--------------------------------------------------------------------------
% Problem 1(a): Gauss-Chebyshev Quadrature 
%--------------------------------------------------------------------------
%%
% Limits of integration 
a = -1.96; 
b = 1.96; 

% Number of evaluation points 
N = [10,20,30,50,75,100,1000,10000,50000]; 

% Standard normal pdf
phi = @(x) normpdf(x);   

% Saving the computed integrals
integral_gcq = zeros(length(N),1);

for k =1:length(N)
    n = N(k);
    % Quadrature nodes
    x = zeros(n, 1); 
    for i=1:n
        val = (2*i-1)/(2*n);
        x(i) = cos(val*pi);
    end 

    % Calculating the integral 
    sum = 0; 
    for i=1:n 
        arg = ((x(i)+1)*(b-a))/2+a; 
        value = phi(arg)*sqrt(1-x(i)^2); 
        sum = sum+value; 
    end 
    integral_gcq(k)= (pi*(b-a))/(2*n)*sum; 
end
%%
%--------------------------------------------------------------------------
% Problem 1(b): Monte Carlo Integration
%--------------------------------------------------------------------------
%%
% Setting seed 
rng(100)

% Limits of integration 
a = -1.96; 
b = 1.96; 

% Number of evaluation points 
n = 10000; 

% Random draws of x_i from Unif(a,b)
x = zeros(n,1);
for i=1:n
    x(i) = unifrnd(a,b);
end 

% Calculating the integral
phi = @(x) normpdf(x);    % standard normal pdf
sum = 0;

for i=1:n
    sum = sum +phi(x(i));
end

integral_mci = (1/n)*(b-a)*sum;
%%
%--------------------------------------------------------------------------
% Problem 1(c): Exact value 
%--------------------------------------------------------------------------
%% 
phi_1 = normcdf(-1.96); 
phi_2 = normcdf(1.96); 
integral_exact = phi_2-phi_1; 

% Differences with the exact result 
diff_gcq = zeros(length(N),1);
for i=1:length(N)
    diff_gcq(i) = abs(integral_exact-integral_gcq(i));
end
diff_mci = abs(integral_exact-integral_mci);