function [theta_mle,beta_mle,Var_beta]=regression_mle(W,n)
% This function compute the maximum likelihood estimation of the model:
% d_{jk}=XB_{jk}+E_{jk}
% where E(E_{jk})=0, Var(E_{jk})=\sigma_{jk}^2.
% Note: This is the special case of mixed effect model d_{jk}=XB_{jk}+ZU_{jk}+E_{jk}
% when Z=0.
% Output:
%      theta_mle: the estimated residual variance.
%      beta_mle: the MLE of B.
%      Var_beta: the estimated marginal variance of B_mle.
[beta_mle,theta_mle,XtX_inv]=sweep_simple_regress(W,n);
Var_beta=diag(XtX_inv)*theta_mle'; % The variance of beta_mle
    