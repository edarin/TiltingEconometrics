function [LOGL, SCORE, INFO] = LOGIT_LOGL(delta,D,X,K,sw)

exp_Xdelta  = exp(X*delta);
LOGL        = -(sum(sw.* D .* (X*delta)) - sum(sw .* log(1+exp_Xdelta)));
SCORE       = -(X'*(sw .* (D - (exp_Xdelta ./ (1+exp_Xdelta)))));
INFO        = -((repmat(sw .* (exp_Xdelta ./ (1+exp_Xdelta).^2),1,K) .* X)'*X);