% Shabab Ahmed 
% ECON 525 HW 2-3

%%
%--------------------------------------------------------------------------
% Problem 1 (a)
%--------------------------------------------------------------------------

% Parameter values 

beta = 0.6; 
A = 20; 
alpha = 0.3; 
delta = 0.5; 

% Grid for capital 
K_max = 12;
K_min = 0;
interval = 0.005; 
Kgrid =(K_min:interval:K_max)';
nbk = length(Kgrid);

% Convergence 
crit = 1; 
eps = 10^-5;


tv = zeros(nbk,1); % matrix to save the transformed value function
V = zeros(nbk,1); % initial value function guess
dr = zeros(nbk,1); % decision rule (will contain indices)

% Value function iteration 

itr = 1; 
while crit>eps
    for i=1:nbk
% checking that we input non-negative values to log
        vals = zeros(nbk,1); 
        for j=1:nbk
            util = A*Kgrid(i)^(alpha)-Kgrid(j)+(1-delta)*Kgrid(i);
            if util > 0
                vals(j) = log(util)+beta.*V(j);
            else 
                vals(j) = -Inf;
            end
        end
        [tv(i), dr(i)] = max(vals);
    end

    crit = max(abs(tv-V)); % Compute convergence criterion
    V = tv; % Update the value function
    itr = itr +1 ;
    disp([num2str(itr), 'th iteration ']);
end

% Plotting V(K)
figure(1)
plot(Kgrid, V, 'bo')
grid on
title('Plot of V(K)'); 
xlabel('K') 
ylabel('V(K)')
axis tight

%%
% Problem 1 b (i)

% Parameter values 

beta = 0.6; 
alpha = 0.3;  

% Grid for capital 
K_max = 12;
K_min = 0;
interval = 0.005; 
Kgrid =(K_min:interval:K_max)';
nbk = length(Kgrid);

% Convergence 
crit = 1; 
eps = 10^-5;


tv = zeros(nbk,1); % matrix to save the transformed value function
V = zeros(nbk,1); % initial value function guess
dr = zeros(nbk,1); % decision rule (will contain indices)

% Value function iteration 

itr = 1; 
while crit>eps
    for i=1:nbk
% checking that we input non-negative values to log
        vals = zeros(nbk,1); 
        for j=1:nbk
            util = Kgrid(i)^(alpha)-Kgrid(j);  
            if util > 0
                vals(j) = log(util)+beta.*V(j);
            else 
                vals(j) = -Inf;
            end
        end
        [tv(i), dr(i)] = max(vals);
    end

    crit = max(abs(tv-V)); % Compute convergence criterion
    V = tv; % Update the value function
    itr = itr +1 ;
    disp([num2str(itr), 'th iteration ']);
end

% Plotting V(K)
figure(2)
plot(Kgrid, V, 'ro')
grid on
title('Plot of V(K)'); 
xlabel('K')
ylabel('V(K)')
axis tight

%%
% Problem 1 (b) (iii)


% Compute the analytical value function. 
A=alpha/(1-alpha*beta)
B=(1/(1-beta))*(log(1-beta*alpha)+beta*A*log(beta*alpha))
V_theory=A*log(Kgrid)+B;

figure(3)
plot(Kgrid,V_theory, 'b')
hold on 
plot(Kgrid, V, 'ro')
grid on
xlabel('K')
ylabel('V(K)')
title('Value Function')
legend('Theoretical','Iterated')
axis tight




%%