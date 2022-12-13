function e = expectation(pi,ai)
    u = 0.5772; 
    if ai == 0
        pi_i = zeros(8,1); 
        pi_j = pi; 
    else 
        pi_i = pi; 
        pi_j = zeros(8,1);
    end 
    exp_term = exp(pi_j-pi_i);
    e = u + log(1+exp_term);
end 