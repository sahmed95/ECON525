% Shabab Ahmed 
% ECON 525 HW 7
%%
%--------------------------------------------------------------------------
% Problem 1(c):  
%--------------------------------------------------------------------------
clear; 
% Parameters 
lambda = 2; 
delta = 2; 
m1 = 1; 
m2 = 1.5; 
phi = 1.5; 
beta  = 0.95; 

% Transition matrix 
b_11 = 0.5; 
b_21 = 0.4; 
B = [b_11, 1-b_11; b_21, 1-b_21];

% Possible states 
states = []; 
m = [m1, m2]; 
a = [0,1]; 
for i=1:2
    for j=1:2
        for k=1:2
            states = [states; [a(i), a(j), m(k)]];
        end
    end
end

states = [0,0,1; 1,0,1; 0,1,1,;1,1,1;0,0,1.5;1,0,1.5;0,1,1.5; 1,1,1.5];
           
%%
% Initial guesses of conditional choice probability 

% Choice probability of entry
choice_prob = (1/9:1/9:8/9);
choice_prob = ones(1,8)*(1/2);

% 2 x 2 identity matrix 
I = eye(8,8); 

% Tolerance level 
tol = 1e-5; 
convergence = 1; 
max_iter = 1000;
iter = 0;

%%
F = transitionmatrix(choice_prob, B, states);
coeff= (I-beta*F); 
pi= profit(lambda, delta, phi, choice_prob, states, 1);
e1_1 = expectation(pi, 1); 
e1_0 = expectation(pi, 0); 
sum_1 = choice_prob'.*(pi+e1_1)+(1-choice_prob)'.*(e1_0); 
V = coeff\sum_1;
enter_prob = F;
exit_prob = F; 
ind_enter = find(states(:,1) == 0);
enter_prob(:,ind_enter) = 0; 
ind_exit = find(states(:,1) == 1); 
exit_prob(:,ind_exit) = 0; 

total_enter_prob = sum(enter_prob,2); 
total_exit_prob = sum(exit_prob, 2);

enter_prob_c = enter_prob./total_enter_prob; 
exit_prob_c = exit_prob./total_exit_prob;


pi_1 = pi; 
pi_0 = zeros(8,1);

V_enter = pi_1+e1_1+beta*(enter_prob_c*V); 
V_exit = pi_0 + e1_0 + beta*(exit_prob_c*V);

denom = exp(V_enter)+exp(V_exit);
psi = exp(V_enter)./denom; 
%%
while convergence > tol 
    iter = iter +1;
    if iter < max_iter
        F = transitionmatrix(choice_prob, B, states);
        coeff= (I-beta*F); 
        pi= profit(lambda, delta, phi, choice_prob, states, 1);
        e1_1 = expectation(pi, 1); 
        e1_0 = expectation(pi, 0); 
        sum_1 = choice_prob'.*(pi+e1_1)+(1-choice_prob)'.*(e1_0); 
        V = coeff\sum_1;
        psi = newchoiceprob(pi, V, F, beta, states, e1_1, e1_0);
        diff = choice_prob-psi';
        convergence = max(abs(diff));
        choice_prob = psi'; 
        disp(convergence)
    end
end 
%%