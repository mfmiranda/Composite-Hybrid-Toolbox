function [beta,gamma,alpha]=UpdateBetaNoOrthog_Michelle(beta,Vbeta,PiMat,Wv,model,wavespecsPC,MCMCspecs,psi_sample)

%%%%   
%%%%        * Updates Betas one-at-a-time 
%%%%        * Assumes non-orthogonal design matrix X
%%%%
%%%%    Input:  (p x K) matrix beta: matrix of sampled betas from previous
%%%%                                MCMC iteration.
%%%%            (p x K) matrix Vbeta: matrix of variances of betans
%%%%            (p x K) matrix PiMat: prior probabilities of nonzero
%%%%                                coefficients
%%%%            (p x K) matrix psi_sample: prior variances for nonzero coefficients
%%%%             Wv:cell array of J or K (p x p) matrices XvX: X'(Sig_jk)^(-1)X
%%%%            cell array of J or K (p x 1) vectors Xvd: X'(Sig_jk)d_jk
%%%%            model: contains model information (design matrix)
%%%%            wavespecsPC: contains number of wavelet coefficients per level
%%%%
%%%%    Output: (p x K) matrix beta, containing samples values of fixed
%%%%                                    effects
%%%%            (p x K) matrix gamma, containing indicators for whether
%%%%                                    coefficient is nonzero
%%%%            (p x K) matrix alpha, containing the posterior probability
%%%%                                of nonzero wavelet coefficients
%%%%
%%%%    Functions needed: UpdateBeta(betans,Vbetans,PI,T);


p=model.p;
K=wavespecsPC.K;

Btilde=Vbeta.*Wv.Xvd;
gamma=NaN(p,K);
alpha=NaN(p,K);
Betansi=NaN(p,K);
    
for i=1:p  %#zhu#% compute beta_i conditional on beta_{-i} here.    
    Bi=NaN(1,K);
    for k=1:K
         Bi(k)=Wv.XvX(i,:,k)*beta(:,k); %#zhu#% these are the stuff to conditional on.
    end 
    Betansi(i,:)=Btilde(i,:)+beta(i,:)-Vbeta(i,:).*Bi;
    [beta(i,:),gamma(i,:),alpha(i,:),]=UpdateBeta_Michelle(Betansi(i,:),Vbeta(i,:),PiMat(i,:),MCMCspecs,psi_sample(i,:));
end 


