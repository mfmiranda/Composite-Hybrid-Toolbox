function [newpsi]=UpdatePsi(model,psiparam,W_trans,K,beta)

KT=model.n;
a0=psiparam.a0;
b0=psiparam.b0;
[Wv_trans,~,~]=GetGCP_trans_Michelle(W_trans,model,0,K);

newpsi=zeros(K,1);
for i=1:K
    M=Wv_trans.dvd(:,i)-2*Wv_trans.Xvd(:,i)'*beta(:,i)+ beta(:,i)'*Wv_trans.XvX(:,:,i)*beta(:,i);
   newpsi(i,1)=1/gamrnd(a0+KT/2,2/(M+2*b0));
end