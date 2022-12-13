function psi = newchoiceprob(pi, V, F,beta, states, e1_1, e1_0)
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
end
