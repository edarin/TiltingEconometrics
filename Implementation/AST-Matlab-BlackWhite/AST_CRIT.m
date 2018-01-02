function [CRIT, FOC, SOC] = AST_CRIT(lambda,D,p_W,p_W_index,t_W,M,NQ,sw)

[phi, phi1, phi2] = AST_PHI(lambda,t_W,p_W_index,NQ);    % compute phi and 1st/2nd derivatives

CRIT    = -sum(sw.* (D .* phi - t_W*lambda) .* (p_W ./ NQ));                % IPT criterion
FOC     = -t_W'*(sw.* (D .* phi1 - 1) .* (p_W ./ NQ));                      % gradient
SOC     = -((repmat(sw .* D .* phi2 .* (p_W ./ NQ),1,1 + M) .* t_W)'*t_W);  % hessian