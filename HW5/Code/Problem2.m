% Shabab Ahmed 
% ECON 525 HW 5
%--------------------------------------------------------------------------
% Problem 2: MLE
%--------------------------------------------------------------------------
%%
clear;
data = readtable("data for hw5.csv");

% Generating x and d
x_data = table2array(data(:,2));
d_data = table2array(data(:,3));

% Length of data
n= length(x_data);

% alpha+beta*x
x_theta = @(theta) theta(1)+theta(2)*x_data;  

% Negative log likelihood
neg_log_l = @(theta) -1*calculate_log_likelihood(x_theta(theta), d_data); 

% Minimizing negative log likelihood
theta_0 = [1,1]; 
theta_hat = fminsearch(neg_log_l, theta_0); 
alpha_hat = theta_hat(1); 
beta_hat = theta_hat(2);
x_theta_hat = x_theta(theta_hat); 
%%
%--------------------------------------------------------------------------
% Problem 2(a): First difference formula
%--------------------------------------------------------------------------
h_vals = [1, 1e-1, 1e-2, 1e-3,1e-4, 1e-5, 1e-6];
se_first_vals = [];
for j = 1:length(h_vals)
    h = h_vals(j);
    alpha_h = alpha_hat+h; 
    beta_h = beta_hat+h; 
    theta_alpha_h = [alpha_h, beta_hat];
    theta_beta_h = [alpha_hat, beta_h]; 
    x_alpha = x_theta(theta_alpha_h);
    x_beta = x_theta(theta_beta_h);

    score_first = zeros(n,2);
    for i =1:n
        likelihood = calculate_log_likelihood(x_theta_hat(i),d_data(i));
        deriv_first_alpha=(calculate_log_likelihood(x_alpha(i), d_data(i))-likelihood)/h;
        deriv_first_beta = (calculate_log_likelihood(x_beta(i), d_data(i))-likelihood)/h;
        score_first(i,:) = [deriv_first_alpha, deriv_first_beta];
    end

    % Variance Covariance matrix
    var_covar_first = (score_first'*score_first)\eye(2);

    % Standard errors
    se_first = sqrt(diag(var_covar_first));
    se_first_vals = [se_first_vals, se_first];
end

diff_first = abs(se_first_vals-se); 
norm_diff_first = zeros(length(h_vals),1);
for i =1:length(h_vals)
    norm_diff_first(i)= norm(diff_first(:,i));
end
%%
%--------------------------------------------------------------------------
% Problem 2(b): Central difference formula
%--------------------------------------------------------------------------
se_central_vals = [];
for j=1:length(h_vals)
    h = h_vals(j);
    alpha_plus = alpha_hat+h; 
    beta_plus = beta_hat+h; 
    alpha_minus = alpha_hat-h;
    beta_minus = beta_hat-h;
    theta_alpha_plus = [alpha_plus, beta_hat];
    theta_alpha_minus = [alpha_minus, beta_hat];
    theta_beta_plus = [alpha_hat, beta_plus]; 
    theta_beta_minus = [alpha_hat, beta_minus];
    x_alpha_plus = x_theta(theta_alpha_plus);
    x_alpha_minus = x_theta(theta_alpha_minus);
    x_beta_plus= x_theta(theta_beta_plus);
    x_beta_minus = x_theta(theta_beta_minus);
    
    % Score 
    score_central = zeros(n,2);
    for i =1:n
        likelihood_alpha_plus = calculate_log_likelihood(x_alpha_plus(i),d_data(i));
        likelihood_alpha_minus = calculate_log_likelihood(x_alpha_minus(i),d_data(i));
        likelihood_beta_plus = calculate_log_likelihood(x_beta_plus(i), d_data(i));
        likelihood_beta_minus = calculate_log_likelihood(x_beta_minus(i), d_data(i));
        deriv_central_alpha = (likelihood_alpha_plus-likelihood_alpha_minus)/(2*h);
        deriv_central_beta = (likelihood_beta_plus-likelihood_beta_minus)/(2*h);
        score_central(i,:) = [deriv_central_alpha, deriv_central_beta];
    end
    
    % Variance Covariance matrix
    var_covar_central = (score_central'*score_central)\eye(2);
    
    % Standard errors
    se_central = sqrt(diag(var_covar_central));
    se_central_vals = [se_central_vals, se_central];
end

diff_central = abs(se_central_vals-se); 
norm_diff_centrals= zeros(length(h_vals),1);
for i =1:length(h_vals)
    norm_diff_centrals(i)= norm(diff_central(:,i));
end
%%
%--------------------------------------------------------------------------
% Problem 2(c): Closed form
%--------------------------------------------------------------------------
pdf = normpdf(x_theta_hat);
cdf_1 = normcdf(x_theta_hat);
cdf_2 = 1-cdf_1;
prod_1 = pdf./cdf_1; 
prod_2 = pdf./cdf_2;
deriv_alpha = d_data.*prod_1-(1-d_data).*prod_2;
deriv_beta = x_data.*deriv_alpha;

% Score 
score = [deriv_alpha, deriv_beta];

% Variance Covariance matrix
var_covar_closed = (score'*score)\eye(2);

% Standard errors
se = sqrt(diag(var_covar_closed));
%%

