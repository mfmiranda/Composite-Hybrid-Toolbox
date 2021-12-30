function [W_trans,Wv_trans,Vbetans_trans,X1,D1]=transformMichelle(D0,scale,model,getW,getWv,mvector)
% This function pre-multiply the scaling parameters to D, X get the new
% model. 

mvector2=repmat(mvector',size(D0,2),1);
scale_newpsi2=repmat(scale.newpsi',size(D0,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(scale.alpha',size(D0,1),1);
%SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector2'.*scale_alpha2));
Ng_lam=(scale_newpsi2.^(-0.5)).*(2.^(mvector2'.*scale_alpha2*0.5));
%Ng_lam is the sqrt(SigmaInverse)
%Ng_lam=(10^9)^(-0.5)*ones(410,1);
K=size(D0,2);
D1=Ng_lam.*D0;   
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

if getWv==1      
    
        [Wv_trans,Vbetans_trans,]=GetGCP_trans_Michelle(W_trans,model,0,K);
 
else
    Wv_trans=[];
    Vbetans_trans=[];
end