function [W_trans]=transformMichelle2(D0,scale,model,getW,mvector)
% This function is equivalent to transformMichelle, except that we don't
% have the psi parameter (innovation variance) on the long memory process

scale.newpsi=ones(size(D0,2),1);
mvector2=repmat(mvector',size(D0,2),1);
scale_newpsi2=repmat(scale.newpsi',size(D0,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(scale.alpha',size(D0,1),1);
%SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector2'.*scale_alpha2));
Ng_lam=(scale_newpsi2.^(-0.5)).*(2.^(mvector2'.*scale_alpha2*0.5));
%Ng_lam=theta_mle0^(-0.5)*size(D0,1);
%Ng_lam=(10^9)^(-0.5)*ones(410,1);
K=size(D0,2);
D1=Ng_lam.*D0;   
%Note:Dteste=diag(Ng_lam(:,1))*D0(:,1) == D1):,1);
if getW==1
   W_trans.XtX=NaN(model.p,model.p,K); 
   W_trans.XtD=NaN(model.p,K);   
   W_trans.DtD=D1'*D1;
   
else
   W_trans=[];
end

X1=NaN(model.n,model.p,K);

for j=1:K    
    Xnj=rowprod(Ng_lam(:,j),model.X);    
    X1(:,:,j)=Xnj;
    if getW==1
       W_trans.XtX(:,:,j)=Xnj'*Xnj;        
       W_trans.XtD(:,j)=Xnj'*D1(:,j);
    end
end

