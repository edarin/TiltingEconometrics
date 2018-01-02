
%-------------------------------------------------------------------------%
%- Monte Carlo Design # 1 : p-score smooth & CEF smooth                  -%                
%-------------------------------------------------------------------------%

rand('state',19);
randn('state',9);
silent = 1;

B           = 5000;
N           = 1000;
Q           = 1/2;
NG          = ones(N,1);
sw          = ones(N,1);

c       = 3;
mu_s    = 0; 
sigma_s = 1;
mu_a    = -1/2; 
sigma_a = 1;
mu_Y    = 0;
sigma_Y = sqrt(3.4823);
alpha0  = 0;
alpha1  = 1/2;
alpha2  = 0;   
sigma_X = 1;
gamma0 = mu_Y - alpha0;

% Analytic standard error calculations
se_AST = 1/10;
se_IPW = 0.1007;
se_CHT = 0.0997;

[mu_Ws, sigma2_Ws, beta0, beta1, beta2] = p_score_AST(mu_s, sigma_s, mu_a, sigma_a, c, Q);


delete('AST_MonteCarlo.log');
diary('AST_MonteCarlo.log');
diary on;

disp('_________________________________________________________________________');
disp('Summary of Monte Carlo Design # 1');
disp('_________________________________________________________________________');
Design_Labels_W = ['mu_s       '; 
                   'sigma_s    '; 
                   'mu_a       '; 
                   'sigma_a    '; 
                   'c          ' ];
Design_Labels_Ws =['mu_Ws      '; 
                   'sigma2_Ws  ';];
Design_Labels_XY =['alpha0     '; 
                   'alpha1     '; 
                   'alpha2     '; 
                   'sigma_X    '; 
                   'mu_Y       '; 
                   'sigma_Y    '];   
Design_Labels_ps =['beta0     '; 
                   'beta1     '; 
                   'beta2     ';]; 

disp('Distribution of W in the study and auxiliary populations');
disp('_________________________________________________________________________');
disp([char(Design_Labels_W) char(num2str([mu_s; sigma_s; mu_a; sigma_a; c]))]);    
disp([char(Design_Labels_Ws) char(num2str([mu_Ws; sigma2_Ws]))]);
disp(' ');              
disp('Distribution of X and Y');
disp('_________________________________________________________________________');
disp([char(Design_Labels_XY) char(num2str([alpha0; alpha1; alpha2; sigma_X; mu_Y; sigma_Y]))]);    
disp(' ');
disp('Implied Propensity Score Coefficients');
disp('_________________________________________________________________________');
disp([char(Design_Labels_ps) char(num2str([beta0; beta1; beta2]))]);    
disp(' ');
disp('Target parameter')
disp('_________________________________________________________________________');
disp([char(['gamma0     ']) char(num2str([gamma0]))]);
disp(' ');
diary off;

MCResults     = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

for b=1:B    
    Ws = random('normal',mu_s*ones(2*N,1), sigma_s*ones(2*N,1));
    i = find((Ws>-c).*(Ws<c));
    Ws = Ws(i); Ws = Ws(1:N);
    
    Wa = random('normal',mu_a*ones(2*N,1), sigma_a*ones(2*N,1));
    i = find((Wa>-c).*(Wa<c));
    Wa = Wa(i); Wa = Wa(1:N);

    D     = (random('unif',zeros(N,1),ones(N,1))<=Q);
    h_W   = [D .* Ws + (1-D) .* Wa];
    DY    = [D .* random('normal',mu_Y*ones(N,1), sigma_Y*ones(N,1))];
    mDX   = (1-D) .* (alpha0 + alpha1*(Wa-mu_Ws) + alpha2*((Wa-mu_Ws).^2-sigma2_Ws) + random('normal',zeros(N,1), sigma_X*ones(N,1)));
    t_W   = [ones(N,1) h_W];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults(b,1:4) = [(gamma_AST - gamma0)/se_AST sqrt(VCOV_gamma_AST/N) (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0)/se_IPW sqrt(VCOV_gamma_IPW/N) (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0)/se_IPW sqrt(VCOV_gamma_CHT/N) (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Median Bias      Std Dev       Median StdErr     MSE      Coverage');          
disp('_________________________________________________________________________');
disp('AST');
disp([median(MCResults(:,1)) std(se_AST*MCResults(:,1)) median(MCResults(:,2)) se_AST*sqrt(mean(MCResults(:,1).^2)) mean(MCResults(:,3)) mean(MCResults(:,4))]);
disp('IPW');
disp([median(MCResults_IPW(:,1)) std(se_IPW*MCResults_IPW(:,1)) median(MCResults_IPW(:,2)) se_IPW*sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([median(MCResults_CHT(:,1)) std(se_CHT*MCResults_CHT(:,1)) median(MCResults_CHT(:,2)) se_CHT*sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 2 : p-score rough & CEF smooth                   -%                
%-------------------------------------------------------------------------%

rand('state',2008);
randn('state',2008);

B           = 5000;
N           = 1000;
Q           = 1/2;
NG          = ones(N,1);
sw          = ones(N,1);


c       = 3;
mu_s    = 0; 
sigma_s = 1;
mu_a    = -1/2; 
sigma_a = sqrt(2/3);
mu_Y    = 0;
sigma_Y = sqrt(2.6580);
alpha0  = 0;
alpha1  = 1/2;
alpha2  = 0;   
sigma_X = 1;
gamma0 = mu_Y - alpha0;

% Analytic standard error calculations
se_AST = 0.0937;
se_IPW = 0.0905;
se_CHT = 0.0925;

[mu_Ws, sigma2_Ws, beta0, beta1, beta2] = p_score_AST(mu_s, sigma_s, mu_a, sigma_a, c, Q);

diary on;
disp('_________________________________________________________________________');
disp('Summary of Monte Carlo Design # 2');
disp('_________________________________________________________________________');
Design_Labels_W = ['mu_s       '; 
                   'sigma_s    '; 
                   'mu_a       '; 
                   'sigma_a    '; 
                   'c          ' ];
Design_Labels_Ws =['mu_Ws      '; 
                   'sigma2_Ws  ';];
Design_Labels_XY =['alpha0     '; 
                   'alpha1     '; 
                   'alpha2     '; 
                   'sigma_X    '; 
                   'mu_Y       '; 
                   'sigma_Y    '];   
Design_Labels_ps =['beta0     '; 
                   'beta1     '; 
                   'beta2     ';]; 

disp('Distribution of W in the study and auxiliary populations');
disp('_________________________________________________________________________');
disp([char(Design_Labels_W) char(num2str([mu_s; sigma_s; mu_a; sigma_a; c]))]);    
disp([char(Design_Labels_Ws) char(num2str([mu_Ws; sigma2_Ws]))]);
disp(' ');              
disp('Distribution of X and Y');
disp('_________________________________________________________________________');
disp([char(Design_Labels_XY) char(num2str([alpha0; alpha1; alpha2; sigma_X; mu_Y; sigma_Y]))]);    
disp(' ');
disp('Implied Propensity Score Coefficients');
disp('_________________________________________________________________________');
disp([char(Design_Labels_ps) char(num2str([beta0; beta1; beta2]))]);    
disp(' ');
disp('Target parameter')
disp('_________________________________________________________________________');
disp([char(['gamma0     ']) char(num2str([gamma0]))]);
disp(' ');
diary off;

MCResults = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

for b=1:B   
    
    Ws = random('normal',mu_s*ones(2*N,1), sigma_s*ones(2*N,1));
    i = find((Ws>-c).*(Ws<c));
    Ws = Ws(i); Ws = Ws(1:N);
    
    Wa = random('normal',mu_a*ones(2*N,1), sigma_a*ones(2*N,1));
    i = find((Wa>-c).*(Wa<c));
    Wa = Wa(i); Wa = Wa(1:N);

    D     = (random('unif',zeros(N,1),ones(N,1))<=Q);
    h_W   = [D .* Ws + (1-D) .* Wa];
    DY    = [D .* random('normal',mu_Y*ones, sigma_Y*ones(N,1))];
    mDX   = (1-D) .* (alpha0 + alpha1*(Wa-mu_Ws) + alpha2*((Wa-mu_Ws).^2-sigma2_Ws) + random('normal',zeros(N,1), sigma_X*ones(N,1)));
    t_W   = [ones(N,1) h_W];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults(b,1:4) = [(gamma_AST - gamma0)/se_AST sqrt(VCOV_gamma_AST/N) (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0)/se_IPW sqrt(VCOV_gamma_IPW/N) (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0)/se_IPW sqrt(VCOV_gamma_CHT/N) (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Median Bias      Std Dev       Median StdErr     MSE      Coverage');          
disp('_________________________________________________________________________');
disp('AST');
disp([median(MCResults(:,1)) std(se_AST*MCResults(:,1)) median(MCResults(:,2)) se_AST*sqrt(mean(MCResults(:,1).^2)) mean(MCResults(:,3)) mean(MCResults(:,4))]);
disp('IPW');
disp([median(MCResults_IPW(:,1)) std(se_IPW*MCResults_IPW(:,1)) median(MCResults_IPW(:,2)) se_IPW*sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([median(MCResults_CHT(:,1)) std(se_CHT*MCResults_CHT(:,1)) median(MCResults_CHT(:,2)) se_CHT*sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 3 : p-score smooth & CEF rough                   -%                
%-------------------------------------------------------------------------%

rand('state',2008);
randn('state',2008);

B           = 5000;
N           = 1000;
Q           = 1/2;
NG          = ones(N,1);
sw          = ones(N,1);

c       = 3;
mu_s    = 0; 
sigma_s = 1;
mu_a    = -0.5; 
sigma_a = 1;
mu_Y    = 0;
sigma_Y = sqrt(1.7496);
alpha0  = 0;
alpha1  = 1/2;
alpha2  = -1;   
sigma_X = 1;
gamma0 = mu_Y - alpha0;

% Analytic standard error calculations
se_AST = 0.1076;
se_IPW = 0.1063;
se_CHT = 0.1309;

[mu_Ws, sigma2_Ws, beta0, beta1, beta2] = p_score_AST(mu_s, sigma_s, mu_a, sigma_a, c, Q);

diary on;
disp('_________________________________________________________________________');
disp('Summary of Monte Carlo Design # 3');
disp('_________________________________________________________________________');
Design_Labels_W = ['mu_s       '; 
                   'sigma_s    '; 
                   'mu_a       '; 
                   'sigma_a    '; 
                   'c          ' ];
Design_Labels_Ws =['mu_Ws      '; 
                   'sigma2_Ws  ';];
Design_Labels_XY =['alpha0     '; 
                   'alpha1     '; 
                   'alpha2     '; 
                   'sigma_X    '; 
                   'mu_Y       '; 
                   'sigma_Y    '];   
Design_Labels_ps =['beta0     '; 
                   'beta1     '; 
                   'beta2     ';]; 

disp('Distribution of W in the study and auxiliary populations');
disp('_________________________________________________________________________');
disp([char(Design_Labels_W) char(num2str([mu_s; sigma_s; mu_a; sigma_a; c]))]);    
disp([char(Design_Labels_Ws) char(num2str([mu_Ws; sigma2_Ws]))]);
disp(' ');              
disp('Distribution of X and Y');
disp('_________________________________________________________________________');
disp([char(Design_Labels_XY) char(num2str([alpha0; alpha1; alpha2; sigma_X; mu_Y; sigma_Y]))]);    
disp(' ');
disp('Implied Propensity Score Coefficients');
disp('_________________________________________________________________________');
disp([char(Design_Labels_ps) char(num2str([beta0; beta1; beta2]))]);    
disp(' ');
disp('Target parameter')
disp('_________________________________________________________________________');
disp([char(['gamma0     ']) char(num2str([gamma0]))]);
disp(' ');
diary off;

MCResults = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

for b=1:B    
    
    Ws = random('normal',mu_s*ones(2*N,1), sigma_s*ones(2*N,1));
    i = find((Ws>-c).*(Ws<c));
    Ws = Ws(i); Ws = Ws(1:N);
    
    Wa = random('normal',mu_a*ones(2*N,1), sigma_a*ones(2*N,1));
    i = find((Wa>-c).*(Wa<c));
    Wa = Wa(i); Wa = Wa(1:N);

    D     = (random('unif',zeros(N,1),ones(N,1))<=Q);
    h_W   = [D .* Ws + (1-D) .* Wa];
    DY    = [D .* random('normal',mu_Y*ones(N,1), sigma_Y*ones(N,1))];
    mDX   = (1-D) .* (alpha0 + alpha1*(Wa-mu_Ws) + alpha2*((Wa-mu_Ws).^2-sigma2_Ws) + random('normal',zeros(N,1), sigma_X*ones(N,1)));
    t_W   = [ones(N,1) h_W];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults(b,1:4) = [(gamma_AST - gamma0)/se_AST sqrt(VCOV_gamma_AST/N) (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0)/se_IPW sqrt(VCOV_gamma_IPW/N) (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0)/se_IPW sqrt(VCOV_gamma_CHT/N) (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
        
end

diary on;
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Median Bias      Std Dev       Median StdErr     MSE      Coverage');          
disp('_________________________________________________________________________');
disp('AST');
disp([median(MCResults(:,1)) std(se_AST*MCResults(:,1)) median(MCResults(:,2)) se_AST*sqrt(mean(MCResults(:,1).^2)) mean(MCResults(:,3)) mean(MCResults(:,4))]);
disp('IPW');
disp([median(MCResults_IPW(:,1)) std(se_IPW*MCResults_IPW(:,1)) median(MCResults_IPW(:,2)) se_IPW*sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([median(MCResults_CHT(:,1)) std(se_CHT*MCResults_CHT(:,1)) median(MCResults_CHT(:,2)) se_CHT*sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 4 : p-score rough & CEF smoth                    -%                
%-------------------------------------------------------------------------%

rand('state',2008);
randn('state',2008);

B           = 5000;
N           = 1000;
Q           = 1/2;
NG          = ones(N,1);
sw          = ones(N,1);

c       = 3;
mu_s    = 0; 
sigma_s = 1;
mu_a    = -0.5; 
sigma_a = sqrt(2/3);
mu_Y    = 0;
sigma_Y = sqrt(0.9253);
alpha0  = 0;
alpha1  = 1/2;
alpha2  = -1;   
sigma_X = 1;
gamma0 = mu_Y - alpha0;

% Analytic standard error calculations
se_AST = 0.0941;
se_IPW = 0.0821;
se_CHT = 0.1192;

[mu_Ws, sigma2_Ws, beta0, beta1, beta2] = p_score_AST(mu_s, sigma_s, mu_a, sigma_a, c, Q);

diary on;
disp('_________________________________________________________________________');
disp('Summary of Monte Carlo Design # 4');
disp('_________________________________________________________________________');
Design_Labels_W = ['mu_s       '; 
                   'sigma_s    '; 
                   'mu_a       '; 
                   'sigma_a    '; 
                   'c          ' ];
Design_Labels_Ws =['mu_Ws      '; 
                   'sigma2_Ws  ';];
Design_Labels_XY =['alpha0     '; 
                   'alpha1     '; 
                   'alpha2     '; 
                   'sigma_X    '; 
                   'mu_Y       '; 
                   'sigma_Y    '];   
Design_Labels_ps =['beta0     '; 
                   'beta1     '; 
                   'beta2     ';]; 

disp('Distribution of W in the study and auxiliary populations');
disp('_________________________________________________________________________');
disp([char(Design_Labels_W) char(num2str([mu_s; sigma_s; mu_a; sigma_a; c]))]);    
disp([char(Design_Labels_Ws) char(num2str([mu_Ws; sigma2_Ws]))]);
disp(' ');              
disp('Distribution of X and Y');
disp('_________________________________________________________________________');
disp([char(Design_Labels_XY) char(num2str([alpha0; alpha1; alpha2; sigma_X; mu_Y; sigma_Y]))]);    
disp(' ');
disp('Implied Propensity Score Coefficients');
disp('_________________________________________________________________________');
disp([char(Design_Labels_ps) char(num2str([beta0; beta1; beta2]))]);    
disp(' ');
disp('Target parameter')
disp('_________________________________________________________________________');
disp([char(['gamma0     ']) char(num2str([gamma0]))]);
disp(' ');
diary off;

MCResults = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

for b=1:B   
    
    Ws = random('normal',mu_s*ones(2*N,1), sigma_s*ones(2*N,1));
    i = find((Ws>-c).*(Ws<c));
    Ws = Ws(i); Ws = Ws(1:N);
    
    Wa = random('normal',mu_a*ones(2*N,1), sigma_a*ones(2*N,1));
    i = find((Wa>-c).*(Wa<c));
    Wa = Wa(i); Wa = Wa(1:N);

    D     = (random('unif',zeros(N,1),ones(N,1))<=Q);
    h_W   = [D .* Ws + (1-D) .* Wa];
    DY    = [D .* random('normal',mu_Y*ones(N,1), sigma_Y*ones(N,1))];
    mDX   = (1-D) .* (alpha0 + alpha1*(Wa-mu_Ws) + alpha2*((Wa-mu_Ws).^2-sigma2_Ws) + random('normal',zeros(N,1), sigma_X*ones(N,1)));
    t_W   = [ones(N,1) h_W];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults(b,1:4) = [(gamma_AST - gamma0)/se_AST sqrt(VCOV_gamma_AST/N) (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0)/se_IPW sqrt(VCOV_gamma_IPW/N) (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0)/se_IPW sqrt(VCOV_gamma_CHT/N) (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Median Bias      Std Dev       Median StdErr     MSE      Coverage');          
disp('_________________________________________________________________________');
disp('AST');
disp([median(MCResults(:,1)) std(se_AST*MCResults(:,1)) median(MCResults(:,2)) se_AST*sqrt(mean(MCResults(:,1).^2)) mean(MCResults(:,3)) mean(MCResults(:,4))]);
disp('IPW');
disp([median(MCResults_IPW(:,1)) std(se_IPW*MCResults_IPW(:,1)) median(MCResults_IPW(:,2)) se_IPW*sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([median(MCResults_CHT(:,1)) std(se_CHT*MCResults_CHT(:,1)) median(MCResults_CHT(:,2)) se_CHT*sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 5a : Qin-Zhang Design (a1)                       -%                
%-------------------------------------------------------------------------%

rand('state',1776);
randn('state',1776);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.1;
BETA2 = 0.1;
gamma0 = 2;

for b=1:B   
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D     = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd   = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (a1), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 5b : Qin-Zhang Design (a2)                       -%                
%-------------------------------------------------------------------------%

rand('state',1777);
randn('state',1777);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.2;
BETA2 = 0.2;
gamma0 = 2;

for b=1:B   
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D     = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd   = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (a2), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 5c : Qin-Zhang Design (a3)                       -%                
%-------------------------------------------------------------------------%

rand('state',1778);
randn('state',1778);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.5;
BETA2 = 0.5;
gamma0 = 2;

for b=1:B   
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1 + 2*W2 + random('normal',zeros(N,1),abs(W2));
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D     = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd   = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (a3), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 6a : Qin-Zhang Design (b1)                       -%                
%-------------------------------------------------------------------------%

rand('state',2068);
randn('state',2068);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.1;
BETA2 = 0.1;
gamma0 = 2;

for b=1:B   
    
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));  
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D    = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd  = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (b1), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 6b : Qin-Zhang Design (b2)                       -%                
%-------------------------------------------------------------------------%

rand('state',2069);
randn('state',2069);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.2;
BETA2 = 0.2;
gamma0 = 2;

for b=1:B   
    
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));  
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D    = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd  = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (b2), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;

%-------------------------------------------------------------------------%
%- Monte Carlo Design # 6c : Qin-Zhang Design (b3)                       -%                
%-------------------------------------------------------------------------%

rand('state',2070);
randn('state',2070);

B           = 1000;
N           = 1000;
NG          = ones(N,1);
sw          = ones(N,1);

MCResults_Lin = zeros(B,4);
MCResults_Qrd = zeros(B,4);
MCResults_IPW = zeros(B,3);
MCResults_CHT = zeros(B,3);

BETA0 = 1;
BETA1 = 0.5;
BETA2 = 0.5;
gamma0 = 2;

for b=1:B   
    
    W1  = random('normal',zeros(N,1),ones(N,1));      
    W2  = ones(N,1) + 0.6*W1 + random('normal',zeros(N,1),ones(N,1));
    Y   = 2 + 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));
    X   = 2*W1.^2 - W2 + 3*W2.^2 + random('normal',zeros(N,1),abs(W2));  
    
    p_W  = (1+exp(-BETA0-BETA1*W1-BETA2*W2)).^-1;       
    D    = (random('unif',zeros(N,1),ones(N,1))<=p_W);
    
    h_W   = [W1 W2];
    DY    = D .* Y;
    mDX   = (1-D) .* X;
    t_W   = [ones(N,1) h_W];
    t_W_Qrd  = [ones(N,1) h_W.^2];
    r_W   = t_W;
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W,mDX,DY,NG,sw,silent);
    MCResults_Lin(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_AST, VCOV_gamma_AST, pi_eff, pi_s, pi_a, exitflag] = AST_ATT(D,r_W,t_W_Qrd,mDX,DY,NG,sw,silent);
    MCResults_Qrd(b,1:4) = [(gamma_AST - gamma0) VCOV_gamma_AST/N (gamma_AST - 1.96*sqrt(VCOV_gamma_AST/N)<=gamma0)*(gamma_AST + 1.96*sqrt(VCOV_gamma_AST/N)>=gamma0) exitflag];
    
    [gamma_IPW, VCOV_gamma_IPW] = IPW_ATT(D,h_W,mDX,DY,NG,sw,silent);
    MCResults_IPW(b,1:3) = [(gamma_IPW - gamma0) VCOV_gamma_IPW/N (gamma_IPW - 1.96*sqrt(VCOV_gamma_IPW/N)<=gamma0)*(gamma_IPW + 1.96*sqrt(VCOV_gamma_IPW/N)>=gamma0)];
       
    [gamma_CHT, VCOV_gamma_CHT] = CHT_ATT(D,[ones(N,1) h_W],mDX,DY,NG,sw);
    MCResults_CHT(b,1:3) = [(gamma_CHT - gamma0) VCOV_gamma_CHT/N (gamma_CHT - 1.96*sqrt(VCOV_gamma_CHT/N)<=gamma0)*(gamma_CHT + 1.96*sqrt(VCOV_gamma_CHT/N)>=gamma0)];    
    
end

diary on;
disp('Qin-Zhang Monte Carlo Design (b3), N = 1000');
disp('Sampling Properties of gamma');
disp('_________________________________________________________________________');
disp('Mean Bias      Variance       IQR-Var       Mean V_Hat     MSE        Coverage');          
disp('_________________________________________________________________________');
disp('AST - Linear');
disp([mean(MCResults_Lin(:,1)) var(MCResults_Lin(:,1)) (iqr(MCResults_Lin(:,1))/1.349).^2  mean(MCResults_Lin(:,2)) sqrt(mean(MCResults_Lin(:,1).^2)) mean(MCResults_Lin(:,3)) mean(MCResults_Lin(:,4))]);
disp('AST - Quadratic');
disp([mean(MCResults_Qrd(:,1)) var(MCResults_Qrd(:,1)) (iqr(MCResults_Qrd(:,1))/1.349).^2  mean(MCResults_Qrd(:,2)) sqrt(mean(MCResults_Qrd(:,1).^2)) mean(MCResults_Qrd(:,3)) mean(MCResults_Qrd(:,4))]);
disp('IPW');
disp([mean(MCResults_IPW(:,1)) var(MCResults_IPW(:,1)) (iqr(MCResults_IPW(:,1))/1.349).^2  mean(MCResults_IPW(:,2)) sqrt(mean(MCResults_IPW(:,1).^2)) mean(MCResults_IPW(:,3))]);
disp('CHT');
disp([mean(MCResults_CHT(:,1)) var(MCResults_CHT(:,1)) (iqr(MCResults_CHT(:,1))/1.349).^2  mean(MCResults_CHT(:,2)) sqrt(mean(MCResults_CHT(:,1).^2)) mean(MCResults_CHT(:,3))]);
disp(' ');
diary off;