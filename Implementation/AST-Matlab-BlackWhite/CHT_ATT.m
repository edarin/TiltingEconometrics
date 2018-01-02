function [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,h_W,mDX,DY,NG,sw)

% This function computes the the average treatment effect on the treated
% under exogenous selection using the CEP estimator of Chen, Hong and
% Tarozzi (2004, 2007). Note that because the ATT moment function is additive
% in the unknown parameter (and just identified) the CHT estimator is
% numerically identical to a flexible parametric imputation type estimator (e.g., Rubin, 1977). 
% We estimate the asymptotic variance of the estimated ATT by treating the CHT 
% estimator "as if" it were a sequential GMM estimator (cf., Hirano and Imbens, 2001).

% INPUTS
% D     : N x 1 vector with ith element equal to 1 if ith unit is treated
%         and zero otherwise
% h_W   : N x M matrix composed of functions of pre-treatment unit characteristics
%         (assume that a constant is in the first column of this matrix)    
% DY    : D*Y, with Y the N x 1 vector of outcomes under treatment
% mDX   : (1-D)*X, with X the N x 1 vector of outcomes under control
% NG    : G x 1 vector with gth row equal to the number of units in the gth cluster 
% sw    : N x 1 vector of known sampling weights

% ----------------------------------------------------------------------------------- %
% - STEP 1 : ORGANIZE DATA                                                          - %
% ----------------------------------------------------------------------------------- %

N       = length(D);      % Number of units in sample  
Ns      = sum(D);         % Number of treated units in the sample
Na      = N - Ns;         % Number of control units in the sample
M       = size(h_W,2);    % Dimension of h_W
G       = length(NG);     % Number of clusters in sample
sw      = sw/mean(sw);    % normalize sample weights to have mean one
Q_hat   = mean(sw .* D);  % Marginal probability of selection

% ----------------------------------------------------------------------------------- %
% - STEP 2 : ESTIMATE ATT/GAMMA                                                     - %
% ----------------------------------------------------------------------------------- %

mD_h_W = repmat(1-D,1,M) .* h_W;                      % Regressor matrix for controls
D_h_W = repmat(D,1,M) .* h_W;                         % Regressor matrix for treated units
PI_star_hat = ((repmat(sw,1,M) .* mD_h_W)'*mD_h_W) \ (repmat(sw,1,M) .* mD_h_W)'*mDX; % CEF of control outcome
gamma_CHT   = sum(sw .* DY)/sum(sw .* D)- sum(sw .* (D_h_W*PI_star_hat))/sum(sw .* D); % Estimate of the ATT

% ----------------------------------------------------------------------------------- %
% - STEP 3 : CALCULATE STANDARD ERRORS                                              - %
% ----------------------------------------------------------------------------------- %

% covariance matrix of moments
m1 = [h_W' .* repmat(sw .* (mDX - mD_h_W*PI_star_hat),1,M)'];               % M x N matrix of m_1 moments
m2 = (sw .* D .* (DY - D_h_W*PI_star_hat - gamma_CHT))';                    % 1 x N matrix of m_2 moments
m  = [m1; m2];                                                              % M x N matrix of moments    

% calculate covariance matrix of moment vector taking into account any within-group dependence/clustering
V_m = zeros(M+1,M+1);
for g = 1:1:G     
    % upper & lower bounds for the g-th group
    n1 = (sum(NG(1:g)) - NG(g)) + 1;                 
    n2 = (sum(NG(1:g)) - NG(g)) + NG(g);                
    
    m_g = sum(m(:,n1:n2),2); 
    V_m = V_m + m_g*m_g'/G;    
end

% inverse of jacobian matrix
M1PI    = ((repmat(sw,1,M) .* mD_h_W)'*mD_h_W)/N;
M2PI    = -mean(repmat(sw,1,M) .* D_h_W);
iM      = [inv(M1PI) zeros(M,1); (1/Q_hat)*M2PI*inv(M1PI) -1/Q_hat];

% estimated asymptotic variance for gamma/ATT
VCOV_gamma_CHT = iM*V_m*iM';
VCOV_gamma_CHT = VCOV_gamma_CHT(end,end);

end

