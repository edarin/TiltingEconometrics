%-------------------------------------------------------------------------%
%- This M file and the accompanying Stata files reproduce                -%
%- estimation results presented in the "Auxiliary-to-Study" paper.       -% 
%- These files are provided "as is". I am unable to assist with their    -%
%- interpretation or use. However please do feel free to e-mail me if    -%
%- you find any mistakes at bryan.graham@nyu.edu. This file generates    -%
%- a MATLAB diary file NLSY_AST_Empirical_Example.log.                   -%             
%-------------------------------------------------------------------------%

clear;

%-------------------------------------------------------------------------%
% Load NLSY79 ESTIMATION SAMPLE                                          -%
%-------------------------------------------------------------------------%

% NOTE: Adjust the directory reference below to point to the directory
% where you have placed the NLSY79_AST_Sample.mat replication dataset file and
% the accompanying M files.
cd('C:\Documents and Settings\bsg1\My Documents\BSG_WORK_19W4th\Research\AST_11Spr\EmpiricalApplication\Created_Data');
load NLSY_AST_Sample.mat;

%-------------------------------------------------------------------------%
% ORGANIZE DATA                                                          -%
%-------------------------------------------------------------------------%

D   = black;
N   = length(D);
N1  = sum(D);
NG  = CountGroupMembers(HHID_79);
G   = length(NG);

% normalize NLSY79 sample weights to have mean one
sw = sample_wgts/mean(sample_wgts);

% take inverse norm of AFQT percentiles
norm_AFQT = ICDF('Norm',AFQT_Adj1/100,0,1);

% ------------------------------------------------------------------------%
% MODEL 1: ADJUST FOR AGE AND AFQT DIFFERENCES FLEXIBLY                  -%
% ------------------------------------------------------------------------%

% vector of moments to balance
h_W = [(yearborn==63)           (yearborn==64)...         
       norm_AFQT                norm_AFQT.^2 ...  
       (norm_AFQT<=-2)          (norm_AFQT<=-1.75)           (norm_AFQT<=-1.5)          (norm_AFQT<=-1.25) ...
       (norm_AFQT<=-1)          (norm_AFQT<=-0.75)           (norm_AFQT<=-0.5)          (norm_AFQT<=-0.25) ...
       (norm_AFQT<=0)           (norm_AFQT<=0.25)            (norm_AFQT<=0.5)           (norm_AFQT<=0.75) ...
      ];      

M = size(h_W,2);                       % number of moments to balance

% set up regression model
DY  = black .* LogWage;
mDX = (1 - black) .* LogWage;
t_W   = [ones(N,1) h_W];
r_W   = t_W;

delete('NLSY_AST_Empirical_Example.log');
diary('NLSY_AST_Empirical_Example.log');
diary off;

% ------------------------------------------------------------------------%
% AST point estimates                                                    -%
% ------------------------------------------------------------------------%
   
[gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,0);
se_gamma      = sqrt(diag(VCOV_gamma_AST)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_AST ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)      

diary on;
disp('AST Point estimates and standard errors for "ATT"')
disp([gamma_AST se_gamma t_gamma pv_gamma]);

% inspect the implicit distribution function estimate
disp('Sum of probability weights and 0 10 25 50 75 90 100 percentiles')
disp('Study sample tilt')
i = find(D==1);
sum(pi_s(i))                               % check weights sum to one
prctile(pi_s(i),[0 10 25 50 75 90 100])    % quantiles of weight distribution
disp('Auxiliary sample tilt')
i = find(D==0);
sum(pi_a(i))                               % check weights sum to one
prctile(pi_a(i),[0 10 25 50 75 90 100])    % quantiles of weight distribution

disp('Check of covariate balance')
disp([sum(repmat(pi_s .* D, 1, M) .* h_W)' sum(repmat(pi_a .* (1-D), 1, M) .* h_W)' sum(repmat(pi_eff, 1, M) .* h_W)']); % verify `exact balance'
diary off;

%------------------------------------------------------------------%
%- plot kernel density estimates of calibrated AFQT distributions -%
%------------------------------------------------------------------%

% find a rule-of-thumb bandwidth value
E_AFQT_eff = sum(pi_eff .* norm_AFQT);                               % efficient estimate of mean BLACK AQFT       
STD_AFQT_eff = sqrt(sum(pi_eff .* (norm_AFQT - E_AFQT_eff).^2));     % efficient estimate of standard deviation of BLACK AQFT       
[AQFT_pi_eff] = sortrows([norm_AFQT pi_eff],1);                      % efficient estimate of BLACK AQFT distribution function 
F_AFQT_eff = cumsum(AQFT_pi_eff(:,2));                               
i = find(F_AFQT_eff>=0.25);                                          % efficient estimate of 0.25 quantile of BLACK AQFT distribution
Q25_AFQT_eff = AQFT_pi_eff(min(i),1);
i = find(F_AFQT_eff>=0.75);                                          % efficient estimate of 0.75 quantile of BLACK AQFT distribution
Q75_AFQT_eff = AQFT_pi_eff(min(i),1);                                  
IQR_AFQT_eff = Q75_AFQT_eff-Q25_AFQT_eff;                            % efficient estimate of IQR of BLACK AQFT distribution       
bw = 0.9*min(STD_AFQT_eff,IQR_AFQT_eff/1.34)*N^(-1/5)                % Silverman rule-of-thumb bandwidth

[F_AFQT_eff,afqt_i]     = ksdensity(norm_AFQT,(-3:0.01:3)','weights',pi_eff,'width',0.5*bw);      % `efficient' estimate of BLACK AFQT density
i = find(D==1);
[F_AFQT_s_raw,afqt_i]   = ksdensity(norm_AFQT(i),(-3:0.01:3)','weights',sw(i),'width',0.5*bw*(N1/N)^(-1/5));        % estimate of actual BLACK AFQT density
i = find(D==0);
[F_AFQT_a,afqt_i]       = ksdensity(norm_AFQT(i),(-3:0.01:3)','weights',pi_a(i),'width',0.5*bw*((N-N1)/N)^(-1/5));  % calibrated estimate of counterfactual WHITE AFQT density
[F_AFQT_a_raw,afqt_i]   = ksdensity(norm_AFQT(i),(-3:0.01:3)','weights',sw(i),'width',0.5*bw*((N-N1)/N)^(-1/5));    % estimate of actual WHITE AFQT density

% plot density estimates
figure;
plot(afqt_i,F_AFQT_s_raw,'-',afqt_i,F_AFQT_a_raw,'-.',afqt_i,F_AFQT_a,':');

% Now compute differences in the study sample CDFs of Y and X (Hourly
% Wages)
for w = 500:250:2000;
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,(1 - black) .* (AvgHourlyWages_90to93<=w),black .* (AvgHourlyWages_90to93<=w),NG,sw,1);
    se_gamma      = sqrt(diag(VCOV_gamma_AST)/G);                % Asymptotic std. errors for gamma_hat
    t_gamma       = gamma_AST ./ se_gamma;                       % t-stats
    pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)      

    diary on;
    disp('AST Point estimates and standard errors for "ATT" of CDF Difference')
    disp(w);
    disp([gamma_AST se_gamma t_gamma pv_gamma]);
    diary off;    
end

% ------------------------------------------------------------------------%
% Compute quantiles of counterfactual distributions and bootstrap         %
% pointwise confidence intervals                                          %
% ------------------------------------------------------------------------%

Quantile_Grid = (3:1:97)';        % grid of quantiles to compute
B = 1000;                         % number of bootstrap replications
BoostrapResults = zeros(B,2*length(Quantile_Grid)); % matrix to store bootstrap results in
HHID_79_uniq    = unique(HHID_79);% G x 1 vector of unique household identifiers

% Compute bootstrap for raw and adjusted quantile differences
for b=1:B
    c = randsample(G,G,true);              % sample households at random with replacement
    NG_b = NG(c);                          % respondents in each sampled households
    N_b = sum(NG_b);                       % total number of sampled households
    i_b = zeros(N_b,1);                    % normalize individual indices                        
    
    % find indices corresponding to all respondents in sampled households
    for g = 1:G
        % actual sample indices for gth group in resample
        n1   = (sum(NG(1:c(g))) - NG(c(g))) + 1;                   
        n2   = (sum(NG(1:c(g))) - NG(c(g))) + NG(c(g));
        
        % resample indices for gth group in resample
        n1_b = (sum(NG_b(1:g)) - NG_b(g)) + 1;                 
        n2_b = (sum(NG_b(1:g)) - NG_b(g)) + NG_b(g);
        
        i_b(n1_b:n2_b) = (n1:n2)';
    end
    
    % Compute AST tilts using bth bootstrap sample
    [gamma_AST_b, VCOV_gamma_AST_b, pi_eff_b, pi_s_b, pi_a_b, exitflag] = AST_ATT(D(i_b),r_W(i_b,:),t_W(i_b,:),mDX(i_b),DY(i_b),NG_b,sw(i_b),1);
        
    % Calculate quantiles of counterfactual distributions using bth bootstrap sample
    [Q_Y_s_b Q_X_s_b Q_X_a_b] = CF_Quantiles(D(i_b),DY(i_b) + mDX(i_b),pi_s_b,pi_a_b,Quantile_Grid);

    % Store bootstrapped results
    BoostrapResults(b,:) = [(Q_Y_s_b-Q_X_a_b)' (Q_Y_s_b-Q_X_s_b)'];
end        

% Plot quantile function of black wage distribution, white wage
% distribution under 'black' characteristics distribution, and white
% wage distribution
[Q_Y_s Q_X_s Q_X_a] = CF_Quantiles(D,LogWage,pi_s,pi_a,Quantile_Grid);
figure;
plot(Quantile_Grid,Q_Y_s,'-',Quantile_Grid,Q_X_s,'-.',Quantile_Grid,Q_X_a,':');

% Compute pointwise confidence intervals for raw and adjusted quantile
% differences
Lower95CI = prctile(BoostrapResults,2.5);
Upper95CI = prctile(BoostrapResults,97.5);
figure;
subplot(1,2,1)
plot(Quantile_Grid,Lower95CI(1:length(Quantile_Grid)),'-k',Quantile_Grid,Upper95CI(1:length(Quantile_Grid)),'-k');
subplot(1,2,2)
plot(Quantile_Grid,Lower95CI(length(Quantile_Grid)+1:end),'-k',Quantile_Grid,Upper95CI(length(Quantile_Grid)+1:end),'-k');

% ------------------------------------------------------------------------%
% IPW point estimates                                                    -%
% ------------------------------------------------------------------------%

% compute IPW point estimates and standard errors
[gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,0);
se_gamma      = sqrt(diag(VCOV_gamma_IPW)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_IPW ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)    

diary on;
disp('IPW Point estimates and standard errors for "ATT"')
disp([gamma_IPW se_gamma t_gamma pv_gamma]);
diary off;

% Now compute differences in the study sample CDFs of Y and X (Hourly
% Wages)
for w = 500:250:2000;
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,(1 - black) .* (AvgHourlyWages_90to93<=w),black .* (AvgHourlyWages_90to93<=w),NG,sw,1);
    se_gamma      = sqrt(diag(VCOV_gamma_IPW)/G);                % Asymptotic std. errors for gamma_hat
    t_gamma       = gamma_IPW ./ se_gamma;                       % t-stats
    pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)      

    diary on;
    disp('IPW Point estimates and standard errors for "ATT" of CDF Difference')
    disp(w);
    disp([gamma_IPW se_gamma t_gamma pv_gamma]);
    diary off;    
end

% ------------------------------------------------------------------------%
% CHT point estimates                                                    -%
% ------------------------------------------------------------------------%

% compute PCEP point estimates and standard errors
[gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
se_gamma      = sqrt(diag(VCOV_gamma_CHT)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_CHT ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)    

diary on;
disp('CHT Point estimates and standard errors for "ATT"')
disp([gamma_CHT se_gamma t_gamma pv_gamma]);
diary off;

% Now compute differences in the study sample CDFs of Y and X (Hourly
% Wages)
for w = 500:250:2000;
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],(1 - black) .* (AvgHourlyWages_90to93<=w),black .* (AvgHourlyWages_90to93<=w),NG,sw);
    se_gamma      = sqrt(diag(VCOV_gamma_CHT)/G);                % Asymptotic std. errors for gamma_hat
    t_gamma       = gamma_CHT ./ se_gamma;                       % t-stats
    pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)      

    diary on;
    disp('CHT Point estimates and standard errors for "ATT" of CDF Difference')
    disp(w);
    disp([gamma_CHT se_gamma t_gamma pv_gamma]);
    diary off;    
end

% ------------------------------------------------------------------------%
% MODEL 2: ADJUST FOR AGE AND AFQT DIFFERENCES LINEARLY                  -%
% ------------------------------------------------------------------------%

% vector of moments to balance
h_W = [(yearborn==63)           (yearborn==64)...         
       norm_AFQT ...    
      ];      

M = size(h_W,2);                       % number of moments to balance

% ------------------------------------------------------------------------%
% AST point estimates                                                    -%
% ------------------------------------------------------------------------%
   
[gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,0);
se_gamma      = sqrt(diag(VCOV_gamma_AST)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_AST ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)      

diary on;
disp('AST Point estimates and standard errors for "ATT"')
disp([gamma_AST se_gamma t_gamma pv_gamma]);
diary off;

% ------------------------------------------------------------------------%
% IPW point estimates                                                    -%
% ------------------------------------------------------------------------%

% compute IPW point estimates and standard errors
[gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,0);
se_gamma      = sqrt(diag(VCOV_gamma_IPW)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_IPW ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)    

diary on;
disp('IPW Point estimates and standard errors for "ATT"')
disp([gamma_IPW se_gamma t_gamma pv_gamma]);
diary off;

% ------------------------------------------------------------------------%
% CHT point estimates                                                    -%
% ------------------------------------------------------------------------%

% compute PCEP point estimates and standard errors
[gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
se_gamma      = sqrt(diag(VCOV_gamma_CHT)/G);                % Asymptotic std. errors for gamma_hat
t_gamma       = gamma_CHT ./ se_gamma;                       % t-stats
pv_gamma      = 2*(1-normcdf(abs(t_gamma)));                 % p-values (for H0: gamma=0)    

diary on;
disp('CHT Point estimates and standard errors for "ATT"')
disp([gamma_CHT se_gamma t_gamma pv_gamma]);
diary off;
