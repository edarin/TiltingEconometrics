function [gamma_IPW, VCOV_gamma_IPW, ps_coef_IPW, VCOV_ps_coef_IPW, pi_IPW, p_score_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent)

% This function estimates the average treatment effects on the treated
% using the inverse probability weighting (IPW) estimator described by 
% Hirano, Imbens and Ridder (2003), Imbens (2004) and also Hirano and
% Imbens (2001). We follow the latter paper in terms of implementation
% details.

% INPUTS
% D      : N x 1 vector with ith element equal to 1 if ith unit in merged
%          sample is from the study population and zero if from the auxiliary
%          population (i.e., D is a "treatment" indicator)
% h_W    : h(W), N x M  matrix of functions of always observed covariates
%          (constant not included)
% mDX    : X,    N x 1  vector of observed outcomes for non-treated/auxiliary units 
% DY     : Y,    N x 1  vector of observed outcomes for treated/study units
% NG     : G x 1 vector with gth row equal to the number of units in the gth cluster 
% sw     : N x 1 vector of known sampling weights
% silent : if silent = 1 display less optimization information and use
%          lower tolerance levels

% OUTPUTS
% gamma_IPW         : IPW estimate of gamma (the ATT)
% VCOV_gamma_IPW    : estimated large sample covariance of gamma
% ps_coef_IPW       : IPW propensity score coefficients
% VCOV_ps_coef_IPW  : estimated large sample covariance for unnormalized
%                     p-score coefficients
% pi_IPW            : Na x 1 vector of auxiliary sample IPW weights
% p_score_IPW       : N x 1 vector of fitted propensity scores

% ----------------------------------------------------------------------------------- %
% - STEP 1 : ESTIMATE PROPENSITY SCORE BY SERIES LOGIT                              - %
% ----------------------------------------------------------------------------------- %

[N M]   = size(h_W);
G       = length(NG);                              % Number of clusters in sample
sw      = sw/mean(sw);                             % normalize sample weights to have mean one

Q_hat = mean(sw .* D);                             % Marginal probability of selection
t_W = [ones(N,1) h_W];                             % Propensity score regressor matrix
[delta_IPW VCOV_delta_IPW] = LOGIT(D,h_W,silent,sw); % "Series" logit p-score estimate

% ----------------------------------------------------------------------------------- %
% - STEP 2 : ESTIMATE ATT BY IPW                                                    - %
% ----------------------------------------------------------------------------------- %

e_V = exp(t_W*delta_IPW);                                % IPW 'weight' p_x/(1-p_x)
tot_e_V = mean(sw .* (1-D) .* e_V);                      % Mean of IPW weights
p_W = e_V ./ (1+e_V);                                    % Fitted propensity score                                              
gamma_IPW = mean(sw .* DY/Q_hat) - mean(sw .* (mDX .* e_V))/tot_e_V; % IPW ATT estimate
ps_coef_IPW = delta_IPW;

% recover auxiliary sample IPW probability weights
pi_IPW = ((1 - D) .* e_V)/(N*tot_e_V);
i      = find(D==0);                                     % Indices for auxiliary/control subsample
pi_IPW = pi_IPW(i);

% ----------------------------------------------------------------------------------- %
% - STEP 3 : FORM LARGE SAMPLE VARIANCE-COVARIANCE ESTIMATES                        - %
% ----------------------------------------------------------------------------------- %

% Form moment vector corresponding to full procedure
m1 = [sw .* (D - Q_hat) sw .* ((1-D) .* e_V - tot_e_V)]';                           % 2   x N matrix of m_1 moments
m2 = (t_W .* repmat(sw .* (D - p_W),1,1+M))';                                       % 1+M x N matrix of m_2 moments
m3 = (sw .* (DY/Q_hat - (mDX .* e_V)/tot_e_V - gamma_IPW))';                        % 1   x N matrix of m_3 moments
m  = [m1; m2; m3];                                                                  % 2 + 1 + M + 1 x N matrix of moments    

% calculate covariance matrix of moment vector taking into account any within-group dependence/clustering
V_m = zeros(2+1+M+1,2+1+M+1);
for g = 1:1:G     
    % upper & lower bounds for the g-th group
    n1 = (sum(NG(1:g)) - NG(g)) + 1;                 
    n2 = (sum(NG(1:g)) - NG(g)) + NG(g);                
    
    m_g = sum(m(:,n1:n2),2); 
    V_m = V_m + m_g*m_g'/G;    
end

% Form Jacobian matrix
M1_rho = -eye(2);                                                                   % 2 x 2 
M2_rho = zeros(1+M,2);                                                              % 1+M x 2
M3_rho = [-mean((sw .* DY)/(Q_hat^2)) mean((sw .*(mDX .* e_V))/(tot_e_V^2))];       % 1 x 2

M1_delta = [zeros(1,1+M); mean(repmat(sw .* ((1-D) .* e_V),1,1+M) .* t_W)];         % 2 x 1 + M;
M2_delta = (repmat(sw .* (- e_V ./ (1+e_V).^2),1,1+M) .* t_W)'*t_W/N;               % 1 + M x 1 + M
M3_delta = mean(repmat(- (sw .*(mDX .* e_V))/tot_e_V,1,1+M) .* t_W);                % 1 x 1 + M;

M3_gamma = -1;

M_full = (N/G)*[M1_rho M1_delta zeros(2,1);...
                M2_rho M2_delta zeros(1+M,1);...
                M3_rho M3_delta M3_gamma];
      
% Calculate asymptotic variance of gamma
VCOV_theta_IPW   = inv(M_full)*V_m*inv(M_full)';                              % Covariance matrix for entire 2(1+M) + 1 parameter (K=1)
VCOV_gamma_IPW   = VCOV_theta_IPW(end,end);                                   % Variance of gamma_hat alone
VCOV_ps_coef_IPW = VCOV_theta_IPW((2+1):(2+1+M),(2+1):(2+1+M));               % Variance of unnormalized propensity score coeffients      
p_score_IPW = p_W;