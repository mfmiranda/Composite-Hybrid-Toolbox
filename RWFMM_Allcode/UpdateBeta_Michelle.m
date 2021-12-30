function [beta,gamma,alpha,lin_shrink]=UpdateBeta_Michelle(betans,Vbetans,PiMat,MCMCspecs,psi_sample)

%%%%    UpdateBeta(betans,Vbetans,PI,T): Updates Fixed effects Beta from f(Beta|theta,pi,T,D) 
%%%%    for a Wavelet-based Functional Mixed Model. (assumes orthogonal
%%%%                                                 design matrix X)
%%%%
%%%%    Input:  (p x K) matrix betans: matrix of non-shrinkage estimates
%%%%                                    for fixed effects beta 
%%%%                                (row i computes the non-shrinkage estimator
%%%%                                 conditioning on the values beta from
%%%%                                 previous run of MCMC).
%%%%            (p x K) matrix Vbetans: matrix of variances of betans
%%%%            (p x K) matrix pi: prior probabilities of nonzero
%%%%                                coefficients
%%%%            (p x K) matrix T: prior variances for nonzero coefficients
%%%%
%%%%    Output: (p x K) matrix beta, containing samples values of fixed
%%%%                                    effects
%%%%            (p x K) matrix gamma, containing indicators for whether
%%%%                                    coefficient is nonzero
%%%%            (p x K) matrix alpha, containing the posterior probability
%%%%                                of nonzero wavelet coefficients
%%%%
%%%%    Functions needed: 
%%%%                        
%%%% %#zhu#% Note: this is the function to compute beta_i conditional on
%%%% beta_{-i}, and the values of these are all contained in betans. I
%%%% suspect the the part of computing betans. 
T=psi_sample./Vbetans;
lin_shrink=T./(T+1);
BF=(T+1).^(-.5).*exp(.5*(betans./sqrt(Vbetans)).^2.*lin_shrink); 
PiMat=min(PiMat,1);%T fix some numerical problemas
PiMat=max(PiMat,0);
O=min((PiMat./(1-PiMat)).*BF,MCMCspecs.maxO); 
alpha=O./(O+1); 
%Make sure this probbility is between zero and one
gamma=binornd(1,alpha);
mu=betans.*lin_shrink;  
sigma=sqrt(Vbetans.*lin_shrink); 
beta=gamma.*normrnd(mu,sigma); %#zhu#% here is the posterior sample for beta.
