% Shabab Ahmed 
% ECON 525 HW 6
%%
%--------------------------------------------------------------------------
% Problem 1(a):  
%--------------------------------------------------------------------------
% Grid for x 
x_min = -1; 
x_max = 1; 
interval = 0.001; 
x_grid = [x_min:interval:x_max]';

% True function 
f = @(x) (1+25*x.^2).^(-1); 

% Degree
n=5; 

% Zeros of T_n(x)
x_zeros = zeros(n,1); 
for i=1:n
    val = (2*i-1)/(2*n)*pi; 
    x_zeros(i) = cos(val); 
end

% Alpha_0 
sum = 0; 
for k =1:n
    sum = sum + f(x_zeros(k)); 
end 
alpha_0 = sum/n; 

% Alpha_i 
alpha_i = zeros(n-1, 1); 
for i =1:n-1
    T_i = @(x) cos(i*acos(x)); 
    sum_i = 0; 
    for k=1:n
        sum_i = sum_i+f(x_zeros(k))*T_i(x_zeros(k)); 
    end 
    alpha_i(i) = (2/n)*sum_i; 
end 

alpha_i = [alpha_0; alpha_i];
f_interp =  0;
for i=1:n
    T_i = @(x) cos((i-1)*acos(x)); 
    f_interp = f_interp+alpha_i(i)*T_i(x_grid); 
end 

figure(1)
plot(x_grid,f(x_grid), 'r-')
hold on
plot(x_grid, f_interp, 'b-')
legend('True function', 'Approximation')
xlabel("x")
ylabel("f(x)")

%%
%--------------------------------------------------------------------------
% Problem 1(b):  
%--------------------------------------------------------------------------
% True function 
one_val = @(x) ones(length(x),1);
values_inner = @(x) [-1*one_val(x), 4*(x-0.2)];
f_inner = @(x) max(values_inner(x),[],2); 
values_outer = @(x) [f_inner(x),one_val(x)]; 
f = @(x) min(values_outer(x),[],2);

% Degree
n = 24; 

% Zeros of T_n(x)
x_zeros = zeros(n,1); 
for i=1:n
    val = (2*i-1)/(2*n)*pi; 
    x_zeros(i) = cos(val); 
end

% Alpha_0 
value_f = f(x_zeros); 
alpha_0 = mean(value_f);

% Alpha_i 
alpha_i = zeros(n-1, 1); 
for i =1:n-1
    T_i = @(x) cos(i*acos(x)); 
    sum_i = 0; 
    for k=1:n
        sum_i = sum_i+f(x_zeros(k))*T_i(x_zeros(k)); 
    end 
    alpha_i(i) = (2/n)*sum_i; 
end 

alpha_i = [alpha_0; alpha_i];
f_interp =  0;
for i=1:n
    T_i = @(x) cos((i-1)*acos(x)); 
    f_interp = f_interp+alpha_i(i)*T_i(x_grid); 
end 


figure(2)
plot(x_grid,f(x_grid), 'r-')
hold on
plot(x_grid, f_interp, 'b-')
legend('True function', 'Approximation')
xlabel("x")
ylabel("f(x)")


