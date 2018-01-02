function [phi, phi1, phi2] = AST_PHI(lambda,t_W,p_W_index,NQ)

% This function evaluates the modified phi(v) function for the logit case
% as described in Graham, Pinto and Egel (2011).

% Coefficients on quadratic extrapolation of phi(v) 
c = -(NQ - 1);
b = NQ + (NQ - 1)*log(1/(NQ - 1));
a = -(NQ - 1)*(1 + log(1/(NQ - 1)) + 0.5*(log(1/(NQ - 1)))^2); 
v_star = log(1/(NQ - 1)); 

% Evaluation of phi(v) and derivatives
v          =  p_W_index + t_W*lambda;
phi        =  (v>v_star) .* (v - exp(-v))   + (v<=v_star) .* (a + b*v + 0.5*c*v.^2);
phi1       =  (v>v_star) .* (1 + exp(-v))   + (v<=v_star) .* (b + c*v);
phi2       =  (v>v_star) .* (- exp(-v))     + (v<=v_star) .* c;


