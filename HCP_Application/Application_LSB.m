
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

% Clusters - Local Spatial Basis
ROI_PC=cell(size(Phi_MatrixNew,1),1);
RV_PC=cell(size(Phi_MatrixNew,1),1);
YC_Matrix_PC=cell(size(Phi_MatrixNew,1),1);

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
    YC_Matrix_PC{i}=PCs(:,1:p)*V(:,1:p)';
end
X0=[];
for i=1:size(Phi_MatrixNew,1)
    X0=[X0,ROI_PC{i}];
end

PC_X0=X0;
    
%% Wavelet Transform on time domain
Y_star=PC_X0';
wavespecs0.wavelet='db1';
wavespecs0.ndim=1;
wavespecs0.compress=1;
padding=107;
Y0=[Y_star,repmat(0,size(Y_star,1),padding)]; %Add padding until numer of basis is even
[Dtrans,wavespecs]=DWT_rows(Y0,wavespecs0);
D=Dtrans';
    
XCov = bsxfun(@minus, DesignMatrix, mean(DesignMatrix));
XCov2 = bsxfun(@rdivide, XCov, std(DesignMatrix));
Xt=[XCov2;repmat(0,padding,size(XCov2,2))];%padding with zeros
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
%-----For LSB
 wavespecsPC.K=size(Y_star,1);
   for kk=1:size(RV_PC,1)
        wavespecsPC.Kj(kk)=size(RV_PC{kk},1);
   end
 wavespecsPC.J=size(RV_PC,1);

       
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
BetaMProjec0=posterior_mean;
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

B=size(Results_MCMC,1);

npc_cumu=0;
BetaDataSpace1=cell(298,1); 

BetaMProjec01=squeeze(Contrast_face_places);
ContrastImageFacePlaces=mean(BetaMProjec01);
for i=1:298
    sa=npc_cumu;
    npc_cumu=npc_cumu+size(RV_PC{i},1);
    indexpc=sa+1:npc_cumu;
    aux1=BetaMProjec01(:,indexpc)*RV_PC{i};
    BetaDataSpace1{i}=aux1;
end

X_BetaDataSpace1=[];

for i=1:298     
    X_BetaDataSpace1=[X_BetaDataSpace1, BetaDataSpace1{i}];   
end
tictoc_projection=toc;
tic
[SimbasResult1,upper_CI, lower_CI]=jointband_simbas_Michelle_sdcorrect(X_BetaDataSpace1,0.05);
tictoc_simbas=toc;
SimbasResult_FacePlaces=SimbasResult1;
Howmany=size(find(SimbasResult1<0.05),2);
Contrast_FacePlaces_Image=mean(X_BetaDataSpace1);
%Saving MCMC Results
save(strcat(path2,'Results_MCMC_WMtask_LSB'),'Results_MCMC','Results_MCMC_psi','MCMC_scalenewpsi','RV_PC','p','K','tictoc_mcmc','tictoc_projection','tictoc_variancecomponents','tictoc_simbas')

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

%save_nii(Contrast,strcat(path2,'Contrast_FacesPlaces_LSB'));

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
NClusters125_LSB=length(IndexClusterSize_125);
NVoxels125_LSB=sum(IndexClusterSize_125);
NClusters64_LSB=length(IndexClusterSize_64);
NVoxels64_LSB=sum(IndexClusterSize_64);
NClusters27_LSB=length(IndexClusterSize_27);
NVoxels27_LSB=sum(IndexClusterSize_27);
BiggCl_LSB=IndexClusterSize_125(1);

save(strcat(path2,'Results_Table3_FacePlaces_task_LSB'),'IntWidth','NClusters125_LSB', 'NVoxels125_LSB','NClusters64_LSB','NVoxels64_LSB','NClusters27_LSB','NVoxels27_LSB','BiggCl_LSB')


%% Generate Cluster Masks for Clusters > 64

for ii=1:NClusters64_LSB
TesteImg=zeros(size(C));
Region=find(Aux1==IndexClusterSize_64(ii));
TesteImg(CC.PixelIdxList{Region(1)})=1;
BigCluster=make_nii(TesteImg);
save_nii(BigCluster,strcat(path1,'Results',sep,'Cluster',num2str(ii),'_LSB'))
end

