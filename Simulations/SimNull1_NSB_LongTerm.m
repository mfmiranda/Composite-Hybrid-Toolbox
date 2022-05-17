% This files runs the simulation results presented 
% in Table 1 and Figure 1 

% Define your toolbox pathaddpath(genpath('tools-for-nifti-and-analyse-image'))

path='/Users/michellemiranda/Documents/MATLAB/';
path2=strcat(path,'SimulationT/');

addpath(genpath(strcat(path,'CodeDWT')))
addpath(genpath(strcat(path,'RWFMM_Allcode')))
addpath(genpath(strcat(path,'tools-for-nifti-and-analyse-image')))

load(strcat(path2,'Sim_PhiMatrix'))
load(strcat(path2,'DataNull_longterm_sim1'))
load(strcat(path2,'X_Design'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation 1 - Composite-hybrid basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spatial basis - project data using ROI basis and then PC
tic

% Wavelet Transform on time domain for both Y and X
Y_star=(Y_MatrixC)';
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
tic
[psi_hat,alpha_hat,ConvergenceMatrix]=EstimateScale_LongMemory(D,nite,model,wavespecs);
toc

scale.newpsi=psi_hat(:,nite);
scale.alpha=alpha_hat(:,nite);
save('Scale_parameters_NSB','psi_hat','alpha_hat')
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
wavespecsPC.J=1;
wavespecsPC.Kj(1)=size(Y_star,1);

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

tic
for i=1:B

[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);
[W_trans2]=transformMichelle2(D,scale,model,1,mvector);
scale.newpsi=UpdatePsi(model,psiparam,W_trans2,K,beta);
[beta,gamma,alpha]=UpdateBetaNoOrthog_Michelle(beta,Vbetans_trans,PiMat,Wv_trans,model,wavespecsPC,MCMCspecs,scale.psi); 
MCMC_beta(i,:,:)=beta;
MCMC_scalenewpsi(i,:)=scale.newpsi;
scale.psi=10*Vbetans_trans;
MCMC_scale_psi(i,:,:)=scale.psi;
i
end
tictoc=toc;
% Results MCMC in basis space
Results_MCMC=reshape(MCMC_beta(burnin:thin:B,:,:),801,p*K);
posterior_mean=shiftdim(mean(MCMC_beta(burnin:thin:i-1,:,:),1),1);
posterior_mean_scale=mean(MCMC_scalenewpsi(burnin:thin:i-1,:),1);

% Projecting back
Results_MCMC_reshape=reshape(Results_MCMC,size(Results_MCMC,1),p,K);
Contrast=MCMC_beta(burnin:thin:i-1,3,:)-MCMC_beta(burnin:thin:i-1,2,:);

X_BetaDataSpace1=shiftdim(Contrast,2); %Contains the MCMC samples in data space



%jointbands
[SimbasResult1,upper_CI, lower_CI]=jointband_simbas_Michelle_sdcorrect((X_BetaDataSpace1)',0.05);
Howmany=size(find(SimbasResult1<0.05),2);

I2=load_nii(strcat(path,'Talairach-labels-2mm.nii'));
C=I2.img(32:63,32:63,26:50);
[dimx,dimy,dimz]=size(C);
ssize=[dimx, dimy,dimz];
PosteriorMeanBeta=mean(X_BetaDataSpace1');


BetaArray=reshape(PosteriorMeanBeta',dimx,dimy,dimz);


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
[iv,jv,kv]=ind2sub(size(C), find(SignalContrast==0) );
ssize=[dimx, dimy,dimz];
MaskVoxels_OutsideContrast=sub2ind(ssize,iv,jv,kv);


MSE=mean((PosteriorMeanBeta'-reshape(SignalContrast,25600,1)).^2);

%MSE=sum(sum(sum((BetaArray-SignalContrast).^2)))/25600;
IntervalWidth=mean(upper_CI-lower_CI);
SimbasResultMask=(SimbasResult1<=0.05);
SimbasMask=reshape(SimbasResultMask',dimx,dimy,dimz);

RP=sum(abs(SimbasMask(MaskVoxelsContrast)==abs(SignalContrast(MaskVoxelsContrast))))/size(MaskVoxelsContrast,1);
FP=sum(SimbasMask(MaskVoxels_OutsideContrast))/size(MaskVoxelsContrast,1);


totaltime=toc;
save(strcat(path2,"ResultsSimNull1_NSB"),'Howmany','IntervalWidth','RP','FP','MSE','totaltime')

