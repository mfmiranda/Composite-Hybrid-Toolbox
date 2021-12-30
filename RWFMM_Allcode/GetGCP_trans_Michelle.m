function [Wv_trans,Vbetans_trans,betans_trans]=GetGCP_trans_Michelle(W_trans,model,Update_betans,K)
%%%%%   Compute generalized cross-products matrices
%%%%%           X'Siginv_jk X       X'Siginv_jk Z       X'Siginv_jk d_jk
%%%%%           Z'Siginv_jk X       Z'Siginv_jk Z       Z'Signinv_jk d_jk
%%%%%           d_jk'Siginv_jk X    d_jk'Siginv_jk Z    d_jk'Siginv_jk d_jk
%%%%%   for (j,k)= (1,1), ..., (J,K_J)
%%%%%   
%%%%%   using Sweep method suggested in Wolfinger, Tobias, and Sall, 1994
%%%%%
%
% Input:    
%           wlevels = (J x 1) vector -- specifies how many wavelet coefficients
%                                   per wavelet level.
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%           theta = (H+c x K) matrix of within-curve variance component
%                           estimates.
%                           
%           W = structure containing:
%
%                 XtX=X'X
%                 XtZ=X'Z
%                 XtD=X'D
%                 ZtD=Z'D
%
%           Wv = structure containing:
% Output:   XvX = cell array with K elements, each containing X' Sigma_{jk}^(-1) X
%           Xvd, dvd, XvZ, ZvZ, Zvd: similarly defined (except make
%               Xvd,Zvd,dvd be matrices of size p x K, m x K, and 1 x K,
%               respectively)
%           L1 = log|Sigma|
%
%  
%
%   Functions needed: repvec.m -- same as Splus function "rep"
%                     sweep.m -- perform sweep operation
%
%   Note: %#zhu#%: In the current version, I removed covQ=covS=0 and 1 cases. 
%   And XvX, XvZ, ZvZ are put in 3-d arrays(rather than cell arrays).
%   Reference: (Sweep operator) Wolfinger, R. 1994, Computing Gaussian
%   Liklihoods and their derivatives for general linear mixed models.

p=model.p;
Xrange=1:p; % Define X, Z, and d partitions of "W" matrix discussed in Wolfinger,
% Tobias, and Sall (1994)



XvX=NaN(p,p,K);
Xvd=NaN(p,K);
dvd=NaN(1,K);


   Drange=(Xrange(end)+1);


    XvX=W_trans.XtX;  
    Xvd=W_trans.XtD;    
    dvd=diag(W_trans.DtD)';

 
%%% Now add part to update betans and Vbetans (check if necessary).
if Update_betans==1
    betans_trans=NaN(p,K);
    Vbetans_trans=NaN(p,K);
    for j=1:K
          XvXk=XvX(:,:,j);          
          temp=[XvXk,Xvd(:,j);Xvd(:,j)',dvd(j)]; %#zhu#% for every j or k, either d will change or both XvX and d will change.
          for ii=1:p
             temp=sweep(temp,ii);
          end
          betans_trans(:,j)=temp(1:p,p+1); %#zhu#% the upper right corner is the estimated MLE. This is not computed as the iterative method stated on Page 188, Morris2006.          
          Vbetans_trans(:,j)=(diag(XvXk)).^(-1);  %%% Vbeta should be computed marginally %#zhu#% Note this is diag(X_i'\Sigma_{j,k}^{-1}X_i).^(-1).            
    end
else
    Vbetans_trans=NaN(p,K);
    for k=1:K
        Vbetans_trans(:,k)=diag(XvX(:,:,k)).^(-1);    
    end
    betans_trans=[];
end

Wv_trans.XvX=XvX;
Wv_trans.Xvd=Xvd;
Wv_trans.dvd=dvd;
