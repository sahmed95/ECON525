function log_l = calculate_log_likelihood(x_theta, d_data)
    n = length(d_data); 
    log_l = 0; 
    for i=1:n 
        cdf = normcdf(x_theta(i));
        val = d_data(i)*log(cdf)+(1-d_data(i))*log(1-cdf);
        log_l = log_l+val; 
    end 
end