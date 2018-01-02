function [Q_Y_s Q_X_s Q_X_a] = CF_Quantiles(D,Z,pi_s,pi_a,Quantile_Grid)

% This function computes quantiles of Y and X in the study population using
% the AST distribution function estimates described in Graham, Pinto and
% Egel (2011). It also computes quantiles of the distribution of X in the
% auxiliary population.

% INPUTS:
% D             : D = 1 if a unit is drawn from the study population and zero
%                 otherwise
% Z             : Z = DY + (1-D)X -- corresponds to observed outcome of interest
% pi_s          : AST study sample tilt
% pi_a          : AST auxiliary sample tilt
% Quantile_Grid : Grid of values to compute corresponding quantiles

% OUTPUTS:
% Q_Y_s         : Quantiles of Y in the study population
% Q_X_s         : Quantiles of X in the study population
% Q_X_a         : Quantiles of X in the auxiliary population

N = length(D);
N1  = sum(D);
numQ = length(Quantile_Grid);

% Compute quantiles of Y distribution in study population
i = find(D==1);
[Y_pi_s] = sortrows([Z(i) pi_s(i)],1);
F_Y_s = cumsum(Y_pi_s(:,2));
Q_Y_s = zeros(numQ,1);

% Compute quantiles of X distribution in study population
i = find(D==0);
[X_pi_s] = sortrows([Z(i) pi_a(i)],1);
F_X_s = cumsum(X_pi_s(:,2));
Q_X_s = zeros(numQ,1);

% Compute quantiles of X distribution in auxiliary population
[X_pi_a] = sortrows([Z(i) ones(N-N1,1)/(N-N1)],1);
F_X_a = cumsum(X_pi_a(:,2));
Q_X_a = zeros(numQ,1);

for j = 1:numQ
    i = find(F_Y_s>=Quantile_Grid(j)/100);
    Q_Y_s(j) = Y_pi_s(min(i),1);    
    i = find(F_X_s>=Quantile_Grid(j)/100);
    Q_X_s(j) = X_pi_s(min(i),1);
    i = find(F_X_a>=Quantile_Grid(j)/100);
    Q_X_a(j) = X_pi_a(min(i),1);    
end

