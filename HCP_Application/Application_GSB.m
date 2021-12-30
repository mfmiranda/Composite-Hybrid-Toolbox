
load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'HCP_Application',sep);
path2=strcat(path1,'Results',sep);
%path3=strcat(path1,'Datasets',sep);

addpath(genpath(strcat(path,sep,'CodeDWT')))
addpath(genpath(strcat(path,sep,'RWFMM_Allcode')))
addpath(genpath(strcat(path,sep,'tools-for-nifti-and-analyse-image')))

load(strcat(path1,'DesignMatrix_100307.mat'))
load(strcat(path1,'InverseTrans_ROI_part_Phi.mat'))
load(strcat(path1,'ExampleSubjectData.mat'))

% Spatial basis - fPCA - global spatial basis

Auxiliar=sum(Y_MatrixC,1);
MaskVoxels_spatial=find(Auxiliar~=0);
Y_MatrixNonZeros=Y_MatrixC(:,MaskVoxels_spatial);

[U0,S0,V0] =svd(Y_MatrixNonZeros,'econ');
VarV0=cumsum(diag(S0).^2)/sum(diag(S0).^2);
idx0=1:1:size(VarV0,1);
PC_X0=Y_MatrixNonZeros*V0;
    
%% Wavelet Transform on time domain
Y_star=PC_X0';
wavespecs0.wavelet='db1';
wavespecs0.ndim=1;
wavespecs0.compress=1;
padding=107;
Y0=[Y_star,repmat(0,size(Y_star,1),padding)]; 
%Add padding until numer of basis is even. This step makes sure that Y=XB+E is the same as
%WY=WXB+WE

[Dtrans,wavespecs]=DWT_rows(Y0,wavespecs0);
D=Dtrans';
    
XCov = bsxfun(@minus, DesignMatrix, mean(DesignMatrix));
XCov2 = bsxfun(@rdivide, XCov, std(DesignMatrix));
Xt=[XCov2;repmat(0,padding,size(XCov2,2))];

[DX,wavespecsX]=DWT_rows(Xt',wavespecs);
model.X=[ones(size(DX,2),1),DX'];
 %% Model set up
[model.n model.p]=size(model.X);
X=model.X;
mvector=[];
    for i=1:wavespecs.J
        aux0=i*ones(wavespecs.Kj(i),1);
        mvector=[mvector;aux0];
    end
W0=GetW_Michelle(model,D);    
[theta_mle0,beta_mle0,Var_beta0]=regression_mle(W0,model.n);

% Estimate scale parameters
tic
[psi_hat,alpha_hat,ConvergenceMatrix]=EstimateScale_LongMemory(D,50,model,wavespecs);
tictoc_variancecomponents=toc;
scale.newpsi=psi_hat(:,50);
scale.alpha=alpha_hat(:,50);
scale.psi=(10^10)*ones(model.p,size(beta_mle0,2));
mvector2=repmat(mvector',size(D,2),1);

% Initial values
Beta_gls=zeros(size(beta_mle0));
for k=1:size(beta_mle0,2)
scale_newpsi2=repmat(scale.newpsi(k,1),size(D,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(scale.alpha(k,1),size(D,1),1);
SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector.*scale_alpha2));
Beta_gls(:,k)=inv(X'*diag(SigmaInv)*X)*X'*diag(SigmaInv)*D(:,k);
end
[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);

% MCMC specifications and hyperparameters
MCMCspecs.nu2_psi_idx=1;
MCMCspecs.minVC=10^(-3);
max_nu2=1e-2; 
extra_scale=1e2;
B=6000;
[p,K]=size(beta_mle0);
MCMC_beta=NaN(B,p,K);
MCMC_scale_psi=NaN(B,p,K);
MCMC_scalenewpsi=NaN(B,K);
   
beta=Beta_gls;
MCMCspecs.maxO=10^20;
gamma=ones(size(beta_mle0));

%Grouping spatial basis
wavespecsPC.K=size(Y_star,1);  
V0Criteria=log2(diag(S0));
wavespecsPC.J=5;
wavespecsPC.Kj(1)=sum(V0Criteria>=18);
wavespecsPC.Kj(2)=sum(V0Criteria<18 & V0Criteria>=17.5);
wavespecsPC.Kj(3)=sum(V0Criteria<17.5 & V0Criteria>=17.2);
wavespecsPC.Kj(4)=sum(V0Criteria<17.2 & V0Criteria>=17);
wavespecsPC.Kj(5)=wavespecsPC.K-sum(wavespecsPC.Kj);
       
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
end
tictoc_mcmc=toc;

burnin=2000;
thin=5;
Results_MCMC=reshape(MCMC_beta(burnin:thin:end,:,:),801,p*K);
Results_MCMC_psi=reshape(MCMC_scale_psi(burnin:thin:end,:,:),801,p*K);
posterior_mean=shiftdim(mean(MCMC_beta(burnin:thin:end,:,:),1),1);
posterior_mean_scale=mean(MCMC_scalenewpsi(burnin:thin:end,:),1);


% Projecting back
tic
Results_MCMC_reshape=reshape(Results_MCMC,size(Results_MCMC,1),p,K);

posterior_mean=squeeze(mean(Results_MCMC_reshape,1));
BetaMProjec0=posterior_mean*V0';
index_2bk=[1;3;5;7];
index_0bk=[9;11;13;15];
index_face=[3,11];
index_tools=[7,15];
index_body=[1,9];
index_places=[5,13];

[p,K]=size(posterior_mean);

%Working memory Contrasts
Contrast_2back_vs_0back_Faces=Results_MCMC_reshape(:,4,:)-Results_MCMC_reshape(:,12,:);
Contrast_2back_vs_0back=mean(Results_MCMC_reshape(:,index_2bk+1,:),2)-mean(Results_MCMC_reshape(:,index_0bk+1,:),2);
Contrast_2back=mean(Results_MCMC_reshape(:,index_2bk+1,:),2);
Contrast_0back=mean(Results_MCMC_reshape(:,index_0bk+1,:),2);

%Category Contrasts
Contrast_body=mean(Results_MCMC_reshape(:,index_body+1,:),2);
Contrast_tools=mean(Results_MCMC_reshape(:,index_tools+1,:),2);
Contrast_face=mean(Results_MCMC_reshape(:,index_face+1,:),2);
Contrast_places=mean(Results_MCMC_reshape(:,index_places+1,:),2);

%Other contrasts
Contrast_face_tools=mean(Results_MCMC_reshape(:,index_face+1,:),2)-mean(Results_MCMC_reshape(:,index_tools+1,:),2);
Contrast_face_places=mean(Results_MCMC_reshape(:,index_face+1,:),2)-mean(Results_MCMC_reshape(:,index_places+1,:),2);
Contrast_body_tools=mean(Results_MCMC_reshape(:,index_body+1,:),2)-mean(Results_MCMC_reshape(:,index_tools+1,:),2);

ProjectV0=V0';
B=size(Results_MCMC,1);

npc_cumu=0;
BetaDataSpace1=cell(298,1); 

BetaMProjec01=squeeze(Contrast_face_places)*ProjectV0;
X_BetaDataSpace1=BetaMProjec01;

tictoc_projection=toc;
tic
[SimbasResult1,upper_CI, lower_CI]=jointband_simbas_Michelle_sdcorrect(X_BetaDataSpace1,0.05);
tictoc_simbas=toc;
%Saving MCMC Results
save(strcat(path2,'Results_MCMC_WMtask_GSB'),'Results_MCMC','Results_MCMC_psi','MCMC_scalenewpsi','V0','p','K','tictoc_mcmc','tictoc_simbas','tictoc_projection','tictoc_variancecomponents')

SimbasResult_FacePlaces=SimbasResult1;
Howmany=size(find(SimbasResult1<0.05),2);
Contrast_FacePlaces_Image=mean(X_BetaDataSpace1);
%save(strcat(path2,'Results_Simbas_FacePlaces_task_CHSB'),'SimbasResult_FacePlaces','upper_CI', 'lower_CI')

% Generate images

I2=load_nii(strcat(path,sep,'Talairach-labels-2mm.nii'));
SimbasResult=SimbasResult_FacePlaces;
C=I2.img;
[dimx, dimy,dimz]=size(C);
[valD numD] = howmany(C);
index_lesspart=find(numD>125);

Size_ROI=numD(index_lesspart(2:end));
AuxSize=cumsum(Size_ROI);
IndexPart=valD(index_lesspart);
IndexPart(1)=[];
SimbasArray=zeros(size(C));
SimbasMask=zeros(size(C));
ContrastImage=zeros(size(C));
SimbasResultMask=(SimbasResult<=0.05);
roi_i=1;
[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
SimbasArray(MaskVoxels)=SimbasResult(1,1:AuxSize(roi_i));
SimbasMask(MaskVoxels)=SimbasResultMask(1,1:AuxSize(roi_i));
ContrastImage(MaskVoxels)=Contrast_FacePlaces_Image(1,1:AuxSize(roi_i));
for roi_i=2:size(IndexPart,1)

[iv,jv,kv]=ind2sub(size(C), find(C==IndexPart(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
SimbasArray(MaskVoxels)=SimbasResult(1,AuxSize(roi_i-1)+1:AuxSize(roi_i));
SimbasMask(MaskVoxels)=SimbasResultMask(1,AuxSize(roi_i-1)+1:AuxSize(roi_i));
ContrastImage(MaskVoxels)=Contrast_FacePlaces_Image(1,AuxSize(roi_i-1)+1:AuxSize(roi_i));

end

SimbasImage1=make_nii(SimbasArray);
SimbasImage2=make_nii(SimbasMask);
Contrast=make_nii(ContrastImage);

save_nii(Contrast,strcat(path2,'Contrast_FacesPlaces_GSB'));

%Clustering
C=SimbasImage2.img;
[iv,jv,kv]=ind2sub(size(C), find(C==1) );
CC = bwconncomp(C,6);
Aux1=zeros(size(CC.PixelIdxList,2),1);
for i=1:size(CC.PixelIdxList,2)
Aux1(i,1)=size(CC.PixelIdxList{i},1);
end
ClusterSizes=sort(Aux1,'descend');
IndexClusterSize_125=ClusterSizes(ClusterSizes>125);
IndexClusterSize_64=ClusterSizes(ClusterSizes>64);
IndexClusterSize_27=ClusterSizes(ClusterSizes>27);

%% Results Table 3
IntWidth=mean(upper_CI-lower_CI);
NClusters125_CHSB=length(IndexClusterSize_125);
NVoxels125_CHSB=sum(IndexClusterSize_125);
NClusters64_CHSB=length(IndexClusterSize_64);
NVoxels64_CHSB=sum(IndexClusterSize_64);
NClusters27_CHSB=length(IndexClusterSize_27);
NVoxels27_CHSB=sum(IndexClusterSize_27);
BiggCl_CHSB=ClusterSizes(1);

save(strcat(path2,'Results_Table3_FacePlaces_task_GSB'),'IntWidth','NClusters125_CHSB', 'NVoxels125_CHSB','NClusters64_CHSB','NVoxels64_CHSB','NClusters27_CHSB','NVoxels27_CHSB','BiggCl_CHSB')


%% Generate Cluster Masks for Clusters > 64

for ii=1:NClusters64_CHSB
TesteImg=zeros(size(C));
Region=find(Aux1==IndexClusterSize_64(ii));
TesteImg(CC.PixelIdxList{Region(1)})=1;
BigCluster=make_nii(TesteImg);
save_nii(BigCluster,strcat(path1,'Results',sep,'Cluster',num2str(ii),'_GSB'))
end


