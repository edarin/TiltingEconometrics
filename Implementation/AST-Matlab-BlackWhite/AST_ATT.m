function [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent)

% This function estimates the average treatment effects on the treated
% using the "auxiliary-to-study tilting" (AST) method described by 
% Graham, Pinto and Egel (2011). The notation mirrors that in the paper
% where possible.

% INPUTS
% D      : N x 1 vector with ith element equal to 1 if ith unit in the merged
%          sample is from the study population and zero if from the auxiliary
%          population (i.e., D is the "treatment" indicator)
% r_W    : r(W), N x 1+L matrix of functions of always observed covariates
%          (constant included) -- these are the propensity score functions
% t_W    : t(W), N x 1+M matrix of functions of always observed covariates
%          (constant included) -- these are the balancing functions
% mDX    : (1-D)*X, N x 1  vector of observed outcomes for non-treated/auxiliary units 
% DY     : D*Y, N x 1  vector of observed outcomes for treated/study units
% NG     : G x 1 vector with gth row equal to the number of units in the gth cluster 
% sw     : N x 1 vector of known sampling weights
% silent : if silent = 1 display less optimization information and use
%          lower tolerance levels

% OUTPUTS
% gamma_AST         : AST estimate of gamma (the ATT)
% VCOV_gamma_AST    : estimated large sample covariance of gamma
% pi_eff            : Semiparametrically efficient estimate of F_s(W) 
% pi_s              : Study sample tilt
% pi_a              : Auxiliary sample tilt 
% exitflag          : equal to 1 if both p-score and tilting optimizations successful  

% Functions called  : LOGIT()              (...LOGIT_LOGL()...)
%                     AST_CRIT()           (...AST_PHI()...)

% ----------------------------------------------------------------------------------- %
% - STEP 1 : ORGANIZE DATA                                                          - %
% ----------------------------------------------------------------------------------- %

N       = length(D);      % Number of units in sample  
Ns      = sum(D);         % Number of units study units in the sample (treated units) 
Na      = N-Ns;           % Number of units auxiliary units in the sample (control units) 
M       = size(t_W,2)-1;  % Dimension of t_W (excluding constant)
L       = size(r_W,2)-1;  % Dimension of r_W (excluding constant)
G       = length(NG);     % Number of clusters in sample
sw      = sw/mean(sw);    % normalize sample weights to have mean one

% ----------------------------------------------------------------------------------- %
% - STEP 2 : ESTIMATE PROPENSITY SCORE PARAMETER BY CMLE                            - %
% ----------------------------------------------------------------------------------- %

[delta_ML VCOV_delta_ML] = LOGIT(D,r_W(:,2:end),silent,sw); % CMLE of p-score coefficients
p_W                      = (1+exp(-r_W*delta_ML)).^-1; % Fitted p-score values
p_W_index                = r_W*delta_ML;               % Fitted p-score index 
NQ                       = sum(sw .* p_W);             % Sum of fitted p-scores
pi_eff                   = (sw .* p_W) ./ NQ;          % Efficient estimate of F_s(W)

% ----------------------------------------------------------------------------------- %
% - STEP 3 : SOLVE FOR AST TILTING PARAMETERS                                       - %
% ----------------------------------------------------------------------------------- %

% Set optimization parameters
if silent == 1
  % Use large scale method with with starting vector of zeros 
  options_delta = optimset('LargeScale','on','GradObj','on','Hessian','on',...
                           'Display','off','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',1000,'MaxIter',1000); 
else
  % Use large scale method with starting vector of zeros and high levels
  % of tolerance; show iteration output
  options_delta = optimset('LargeScale','on','GradObj','on','Hessian','on',...
                           'Display','iter','TolFun',1e-25,'TolX',1e-12,'MaxFunEvals',10000,'MaxIter',10000); 
end

lambda_SV = zeros(1+M,1); % use vector of zeros as starting values
  
% Compute lamba_s_hat (study or treated sample tilting parameters)
f_s_ast_crit = @(x)AST_CRIT(x,D,p_W,p_W_index,t_W,M,NQ,sw);  
[lambda_s_hat, fval, exitflag1] = fminunc(f_s_ast_crit, lambda_SV, options_delta); 
p_W_s = (1+exp(-r_W*delta_ML-t_W*lambda_s_hat)).^-1; % study sample tilted p-score
pi_s = D .* pi_eff ./ p_W_s;                         % study sample tilt 

% Compute lamba_a_hat (auxiliary or non-treated sample tilting parameters)
f_a_ast_crit = @(x)AST_CRIT(x,1-D,p_W,-p_W_index,t_W,M,NQ,sw);  
[lambda_a_hat, fval, exitflag2] = fminunc(f_a_ast_crit, lambda_SV, options_delta); 
lambda_a_hat = -lambda_a_hat;
p_W_a = (1+exp(-r_W*delta_ML-t_W*lambda_a_hat)).^-1; % auxiliary sample tilted p-score
pi_a = (1-D) .* pi_eff ./ (1-p_W_a);                 % auxiliary sample tilt 
exitflag = (min(exitflag1,exitflag2)>0);             % equal one if both optimizations successful, zero otherwise

% ----------------------------------------------------------------------------------- %
% - STEP 4 : SOLVE FOR AST ESTIMATE OF GAMMA (i.e., ATT)                            - %
% ----------------------------------------------------------------------------------- %

% AST estimate of gamma -- the ATT %
gamma_AST = sum(sw .* p_W .* ((D ./ p_W_s) .* DY - (1-D) ./ (1-p_W_a) .* mDX))/NQ;

% ----------------------------------------------------------------------------------- %
% - STEP 5 : FORM LARGE SAMPLE VARIANCE-COVARIANCE ESTIMATES                        - %
% ----------------------------------------------------------------------------------- %

% Form moment vector corresponding to full three step procedure
m1 = (repmat(sw .* (D - p_W),1,1+L) .* r_W)';                                        % 1+L x N matrix of m_1 moments (logit scores)
m2 = (repmat(sw .* ((1 - D) ./ (1 - p_W_a) - 1) .* p_W,1,1+M) .* t_W)';              % 1+M x N matrix of m_2 moments
m3 = (repmat(sw .* (D ./ p_W_s - 1) .* p_W,1,1+M) .* t_W)';                          % 1+M x N matrix of m_3 moments
m4 = (sw .* p_W .* ((D ./ p_W_s) .* DY - ((1-D) ./ (1-p_W_a)) .* (mDX+gamma_AST)))'; % 1 x N matrix of m_4 moments
m  = [m1; m2; m3; m4];                                                               % 1 + L + 2(1 + M) + 1 x N matrix of moments    

% calculate covariance matrix of moment vector taking into account any within-group dependence/clustering
V_m = zeros(1+L+2*(1+M)+1,1+L+2*(1+M)+1);
for g = 1:1:G     
    % upper & lower bounds for the g-th group
    n1 = (sum(NG(1:g)) - NG(g)) + 1;                 
    n2 = (sum(NG(1:g)) - NG(g)) + NG(g);                
    
    m_g = sum(m(:,n1:n2),2); 
    V_m = V_m + m_g*m_g'/G;    
end

% Form Jacobian matrix for entire parameter: theta = (rho, delta, gamma)
e_V  = exp(r_W*delta_ML);
e_Va = exp(r_W*delta_ML+t_W*lambda_a_hat);
e_Vs = exp(r_W*delta_ML+t_W*lambda_s_hat);

M1_delta = (repmat(sw .* (- e_V ./ (1 + e_V).^2),1,1+L) .* r_W)'*r_W/N;                                               % 1 + L x 1 + L
M2_delta = (repmat(sw .* ((1 - D) ./ (1 - p_W_a) - 1) .* (e_V ./ (1 + e_V).^2),1,1+M) .* t_W)'*r_W/N;                 % 1 + M x 1 + L     
M3_delta = (repmat(sw .* (D ./ p_W_s - 1) .* (e_V ./ (1 + e_V).^2),1,1+M) .* t_W)'*r_W/N;                             % 1 + M x 1 + L     
M4_delta = (sw .* (e_V ./ (1 + e_V).^2) .* ((D ./ p_W_s) .* DY - ((1 - D) ./ (1 - p_W_a)) .* (mDX + gamma_AST)))'*r_W/N;% 1     x 1 + L    

M2_lambda_a = (repmat(sw .* ((1 - D) ./ (1 - p_W_a).^2) .* p_W .* (e_Va ./ (1 + e_Va).^2),1,1+M) .* t_W)'*t_W/N;      % 1 + M x 1 + M
M4_lambda_a = (      -sw .* ((1 - D) ./ (1 - p_W_a).^2) .* p_W .* (mDX+gamma_AST) .* (e_Va ./ (1 + e_Va).^2))' * t_W/N;% 1     x 1 + M                                    

M3_lambda_s = (repmat(-sw .* (D ./ p_W_s.^2) .* p_W .* (e_Vs ./ (1 + e_Vs).^2),1,1+M) .* t_W)'*t_W/N;                 % 1 + M x 1 + M
M4_lambda_s = (       -sw .* (D ./ p_W_s.^2) .* p_W .* DY .* (e_Vs ./ (1 + e_Vs).^2))' * t_W/N;                       % 1     x 1 + M

M4_gamma = -NQ/N;                                                                                               % 1     x 1  

M_hat = (N/G)*[M1_delta         zeros(1+L,1+M)            zeros(1+L,1+M)            zeros(1+L,1); ...
               M2_delta         M2_lambda_a               zeros(1+M,1+M)            zeros(1+M,1); ...
               M3_delta         zeros(1+M,1+M)            M3_lambda_s               zeros(1+M,1); ...
               M4_delta         M4_lambda_a               M4_lambda_s               M4_gamma];
       
VCOV_gamma_AST  = inv(M_hat)*V_m*inv(M_hat)';
VCOV_gamma_AST  = VCOV_gamma_AST(end,end);       