function ln_L = Likelihood(x_data, i_data, Pi_0, Pi_1, P)
    % Function to calculate the likelihood of the data
    len_data = length(x_data); 
    ln_L= 0;
    for i=2:len_data
        i_t = i_data(i);
        x_t = x_data(i);
        p_it = P(x_t+1)*i_t +(1-i_t)*(1-P(x_t+1));
        if i_data(i-1) == 1
            p_xt = Pi_1(1, x_t+1);
        else 
            p_xt = Pi_0(x_data(i-1)+1,x_t+1);
        end 
        ln_L = ln_L+log(p_xt)+log(p_it);
    end 
end 
