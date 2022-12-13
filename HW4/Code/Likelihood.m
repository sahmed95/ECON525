function L = Likelihood(x_data, i_data, p_hat_0, p_hat_1, V_0, V_1)
    % Function to calculate the likelihood of the data
    P_i = zeros(11,1); 
    for i=1:11
        num = exp(V_1(i));
        denom = exp(V_0(i))+num;
        P_i(i) = num/denom; 
    end 
    len_data = length(x_data); 
    L = 0;
    
    for i=2:(len_data)
        i_t = i_data(i);
        x_t = x_data(i);
        x_t1 = x_data(i-1);
        p_it = P_i(x_t+1)*i_t +(1-i_t)*(1-P_i(x_t+1));
        if i_data(i-1) == 1
            p_xt = p_hat_1(x_t1+1,x_t+1);
        else 
            p_xt = p_hat_0(x_t1+1,x_t+1);
        end
        L = L+log(p_xt)+log(p_it);
    end 
end 
