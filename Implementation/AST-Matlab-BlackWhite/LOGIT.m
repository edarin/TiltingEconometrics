function [delta_ML VCOV_delta_ML] = LOGIT(D,X,silent,sw);

% This function computes the ML estimates of the logit binary choice model:
% Pr(D=1|X=x)= exp(x'delta)/[1+exp(x'delta)]

% INPUTS
% D      : N x 1 vector of binary outcomes
% X      : X is a N x K matrix of covariates (without a constant)
% silent : when silent = 1 optimization output is suppressed and
%          optimization is by Fisher-scoring with lower tolerances.
%          Otherwise optimization starts with a few Fisher-scoring steps
%          and then switches to a quasi-Newton search.         
% sw     : N x 1 vector of known sampling weights (assumed to have mean
%          one)

% OUTPUTS
% gamma_ML         : ML estimates of logit coefficients 
% VCOV_delta_ML    : large sample covariance of estimates (= inverse Infomation)

% Functions called : LOGIT_LOGL()

[N K] = size(X);        % Number of observations and covariates
X = [ones(N,1) X];      % Add a constant to the regressor matrix
K = K + 1;
f_logit_logl = @(x)LOGIT_LOGL(x, D, X, K, sw);  % define objective function

% Set optimization parameters
if silent == 1      
    % Use fisher-scoring with lower tolerances
    options_delta = optimset('LargeScale','on','GradObj','on','Hessian','on',...
                             'Display','off','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',1000,'MaxIter',1000);
    delta_SV = (X'*X) \ (X'*D);                                 % Linear probability starting values for logit          
else
    % Take a few fisher-scoring steps to get starting values, then do a
    % quasi-newton search
    options_delta = optimset('LargeScale','on','GradObj','on','Hessian','on',...
                             'Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',10,'MaxIter',10); 
    
    delta_SV = fminunc(f_logit_logl,  zeros(K,1), options_delta); 
    
    options_delta = optimset('LargeScale','off','GradObj','on','Hessian','off',...
                             'Display','iter','TolFun',1e-25,'TolX',1e-12,'MaxFunEvals',10000,'MaxIter',10000); 
    
end                                      

[delta_ML,LOGL,exitflag,output,SCORE,INFO] = fminunc(f_logit_logl, delta_SV, options_delta);  % estimate delta_ML    
VCOV_delta_ML = inv(INFO);                                 % Use the inverse information to estimate the large sample variance of delta_ML