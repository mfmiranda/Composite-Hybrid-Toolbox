clear

load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'Simulations',sep);
path2=strcat(path,sep,'Simulations',sep,'Replications_100',sep);
path3=strcat(path1,'Datasets',sep);

addpath(genpath(strcat(path,'CodeDWT')))
addpath(genpath(strcat(path,'RWFMM_Allcode')))
addpath(genpath(strcat(path,'tools-for-nifti-and-analyse-image')))


MSEArray=zeros(100,1);
FPArray=zeros(100,1);
RPArray=zeros(100,1);
IntervalWidthArray=zeros(100,1);
load(strcat(path1,'Sim_PhiMatrix'))
load(strcat(path1,'X_Design'))

for dataset=1:100
filename=strcat(path3,'DataLongTermNoise_',num2str(dataset));
load(filename)


% Spatial basis - project data using ROI basis and then PC

ROI_PC=cell(size(Phi_MatrixNew,1),1);
RV_PC=cell(size(Phi_MatrixNew,1),1);

for i=1:size(Phi_MatrixNew,1)
    ind_roi=find(Phi_MatrixNew(i,:)==1);
    X=double(Y_MatrixC(:,ind_roi));
    [U,S,V] =svd(X,'econ');
    VarV=cumsum(diag(S).^2)/sum(diag(S).^2);
    idx=1:1:size(VarV,1);
    p=min(idx(VarV>=0.5));
    PCs=X*V;
    ROI_PC{i}=PCs(:,1:p);
    RV_PC{i}=V(:,1:p)';
end

X0=[];
for i=1:size(Phi_MatrixNew,1)
    X0=[X0,ROI_PC{i}];
end
 %X0 is the matrix of local features

[U0,S0,V0] =svd(X0,'econ');
VarV0=cumsum(diag(S0).^2)/sum(diag(S0).^2);
idx0=1:1:size(VarV0,1);
p0=min(idx0(VarV0>=0.999));
PC_X0=X0*V0(:,1:p0); %PC_X0 is the matrix of global features

% Wavelet Transform on time domain for both Y and X
Y_star=PC_X0';

wavespecs0.wavelet='db1';
wavespecs0.ndim=1;
wavespecs0.compress=1;
Y0=[Y_star,repmat(0,size(Y_star,1),11)]; %Add padding until numer of basis is odd
[Dtrans,wavespecs]=DWT_rows(Y0,wavespecs0);
D=Dtrans'; %Data in the wavelet space
Xt2=[Xt;repmat(0,11,size(Xt,2))];
[DX,wavespecsX]=DWT_rows(Xt2',wavespecs); %Add padding until numer of basis is odd
model.X=[ones(size(DX,2),1),DX']; %Design matrix on the wavelet space
   
% Running the model
[model.n model.p]=size(model.X);
X=model.X;
mvector=[];
    for i=1:wavespecs.J
        aux0=i*ones(wavespecs.Kj(i),1);
        mvector=[mvector;aux0];
    end

W0=GetW_Michelle(model,D);    
[theta_mle0,beta_mle0,Var_beta0]=regression_mle(W0,model.n);

% Estimate long memory and scale parameters
% Done separately
nite=25;
[psi_hat,alpha_hat,ConvergenceMatrix]=EstimateScale_LongMemory(D,nite,model,wavespecs);
scale.newpsi=psi_hat(:,nite);
scale.alpha=alpha_hat(:,nite);
scale.psi=(10^3)*ones(model.p,size(beta_mle0,2));

mvector2=repmat(mvector',size(D,2),1);

Beta_gls=zeros(size(beta_mle0));
for k=1:size(beta_mle0,2)
scale_newpsi2=repmat(psi_hat(k,nite)',size(D,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(alpha_hat(k,nite)',size(D,1),1);
SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector.*scale_alpha2));
Beta_gls(:,k)=inv(X'*diag(SigmaInv)*X)*X'*diag(SigmaInv)*D(:,k);
end
[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);


B=6000; %Number of total MCMC samples
burnin=2000; %Number of burnin samples
thin=5; %MCMC thinning

[p,K]=size(beta_mle0);
MCMC_beta=NaN(B,p,K);
MCMC_scale_psi=NaN(B,p,K);
MCMC_scalenewpsi=NaN(B,K);
    
beta=Beta_gls;
gamma=ones(size(beta_mle0));

MCMCspecs.maxO=10^20;
MCMCspecs.nu2_psi_idx=1;
MCMCspecs.minVC=10^(-3);
max_nu2=1e-2; 
extra_scale=1e2;


%Grouping spatial basis

%-----For CHSB
wavespecsPC.K=size(Y_star,1);
V0Criteria=log2(diag(S0(1:p0,1:p0)));
wavespecsPC.J=4;
wavespecsPC.Kj(1)=sum(V0Criteria>=8.5);
wavespecsPC.Kj(2)=sum(V0Criteria>=7.5 & V0Criteria<8.5);
wavespecsPC.Kj(3)=sum(V0Criteria>=7 & V0Criteria<7.5);
wavespecsPC.Kj(4)=wavespecsPC.K-sum(wavespecsPC.Kj);

 %For L  
       
if MCMCspecs.nu2_psi_idx==1 % 1 indicating that nu2_psi depend on a and j.  
     meanop=uneqkron(wavespecsPC.Kj)*diag(1./wavespecsPC.Kj);
    nu2.psi=min(2./(Var_beta0*meanop)/extra_scale,max_nu2);    
else                     % 0 indicating that nu2_psi depend on j and k.     
    nu2.psi=min(2./mean(Var_beta0)'/extra_scale,max_nu2);
end

PiMat=ones(size(Vbetans_trans));

psiparam.a0=2;
psiparam.b0=2;


for i=1:B

[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);
[W_trans2]=transformMichelle2(D,scale,model,1,mvector);
scale.newpsi=UpdatePsi(model,psiparam,W_trans2,K,beta);
[beta,gamma,alpha]=UpdateBetaNoOrthog_Michelle(beta,Vbetans_trans,PiMat,Wv_trans,model,wavespecsPC,MCMCspecs,scale.psi); 
MCMC_beta(i,:,:)=beta;
MCMC_scalenewpsi(i,:)=scale.newpsi;
scale.psi=10*Vbetans_trans;
MCMC_scale_psi(i,:,:)=scale.psi;
end

% Results MCMC in basis space
Results_MCMC=reshape(MCMC_beta(burnin:thin:end,:,:),801,p*K);
posterior_mean=shiftdim(mean(MCMC_beta(burnin:thin:end,:,:),1),1);
posterior_mean_scale=mean(MCMC_scalenewpsi(burnin:thin:end,:),1);

% Projecting back
Results_MCMC_reshape=reshape(Results_MCMC,size(Results_MCMC,1),p,K);
Contrast=Results_MCMC_reshape(:,3,:)-Results_MCMC_reshape(:,2,:);
npc_cumu=0;
BetaDataSpace1=cell(size(Phi_MatrixNew,1),1); 
ProjectV0=(V0(:,1:p0))';
BetaMProjec01=squeeze(Contrast)*ProjectV0;
for i=1:size(Phi_MatrixNew,1)
    sa=npc_cumu;
    npc_cumu=npc_cumu+size(RV_PC{i},1);
    indexpc=sa+1:npc_cumu;
    aux1=BetaMProjec01(:,indexpc)*RV_PC{i};
    BetaDataSpace1{i}=aux1;
end

X_BetaDataSpace1=[]; %Contains the MCMC samples in data space

for i=1:size(Phi_MatrixNew,1) 
    X_BetaDataSpace1=[X_BetaDataSpace1, BetaDataSpace1{i}];
end

%jointbands
[SimbasResult1,upper_CI, lower_CI]=jointband_simbas_Michelle_sdcorrect(X_BetaDataSpace1,0.05);
Howmany=size(find(SimbasResult1<0.05),2);

I2=load_nii(strcat(path,sep,'Talairach-labels-2mm.nii'));

C=I2.img(32:63,32:63,26:50);
[dimx,dimy,dimz]=size(C);
PosteriorMeanBeta=mean(X_BetaDataSpace1);
BetaArray=zeros(size(C));
roi_i=1;
[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
BetaArray(MaskVoxels)=PosteriorMeanBeta(1,1:AuxSize(roi_i));

for roi_i=2:size(IndexPart,1)
[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
 ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
BetaArray(MaskVoxels)=PosteriorMeanBeta(1,AuxSize(roi_i-1)+1:AuxSize(roi_i));
end

[iv2,jv2,kv2]=ind2sub(size(C), find(C==IndexPart(8)) );
MaskVoxels2=sub2ind(ssize,iv2,jv2,kv2);
Signal2=zeros(dimx,dimy,dimz);
Signal2(MaskVoxels2)=1;

[iv3,jv3,kv3]=ind2sub(size(C), find(C==IndexPart(3)) );
MaskVoxels3=sub2ind(ssize,iv3,jv3,kv3);
Signal3=zeros(dimx,dimy,dimz);
Signal3(MaskVoxels3)=1;

SignalContrast=Signal3-Signal2;
[ivc,jvc,kvc]=ind2sub(size(C), find(C==IndexPart(3)| C==IndexPart(8)));
MaskVoxelsContrast=sub2ind(ssize,ivc,jvc,kvc);
MSE=sum(sum(sum((BetaArray-SignalContrast).^2)))/25600;

%False Positives and false negatives
SimbasArray=zeros(size(C));
SimbasMask=zeros(size(C));
SimbasResultMask=(SimbasResult1<=0.05);
roi_i=1;
[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
SimbasMask(MaskVoxels)=SimbasResultMask(1,1:AuxSize(roi_i));

for roi_i=2:size(IndexPart,1)

[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
SimbasMask(MaskVoxels)=SimbasResultMask(1,AuxSize(roi_i-1)+1:AuxSize(roi_i));

end

[iv,jv,kv]=ind2sub(size(C), find(SignalContrast==0) );
ssize=[dimx, dimy,dimz];
MaskVoxels_OutsideContrast=sub2ind(ssize,iv,jv,kv);
%Detecting true signal - RP=real positives
RP=sum(abs(SimbasMask(MaskVoxelsContrast)==abs(SignalContrast(MaskVoxelsContrast))))/size(MaskVoxelsContrast,1);
FP=sum(SimbasMask(MaskVoxels_OutsideContrast))/size(MaskVoxelsContrast,1);
IntervalWidth=mean(upper_CI-lower_CI);

MSEArray(dataset,1)=MSE;
FPArray(dataset,1)=FP;
RPArray(dataset,1)=RP;
IntervalWidthArray(dataset,1)=IntervalWidth;
clearvars -except MSEArray FPArray RPArray IntervalWidthArray path3 path2 path Xt  Mac sep Phi_MatrixNew IndexPart AuxSize
end
save(strcat(path2,'ResultsSim100_CHSB_LongTerm'),'MSEArray','FPArray', 'RPArray', 'IntervalWidthArray')


