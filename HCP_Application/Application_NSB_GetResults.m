%This function will combine the results of each split file 
%running the non spatial basis model
load(strcat(path1,'ExampleSubjectData.mat'))
Auxiliar=sum(Y_MatrixC,1);
MaskVoxels_spatial=find(Auxiliar~=0);
Y_MatrixNonZeros=Y_MatrixC(:,MaskVoxels_spatial);
Contrast_MCMC_NoBasis=zeros(401,size(Y_MatrixNonZeros,2));
Contrast_face_places=zeros(401,size(Y_MatrixC,2));
for k=1:5
    load(strcat('Results_MCMC_NoBasis',num2str(k)));
    Contrast_MCMC_NoBasis(:,(k-1)*22877+1:k*22877)=Results_MCMC;
end

[SimbasResult1,upper_CI, lower_CI]=jointband_simbas_Michelle_sdcorrect(Contrast_MCMC_NoBasis,0.05);

SimbasResult_FacePlaces=SimbasResult1;
Howmany=size(find(SimbasResult1<0.05),2);
Contrast_FacePlaces_Image=mean(Contrast_MCMC_NoBasis);
%Make Images
Aux0=SimbasResult1<0.05;
SimbasArray=zeros(1,size(Y_MatrixC,2));
SimbasArray(1,MaskVoxels_spatial)=Aux0;
ContrastArray=zeros(1,size(Y_MatrixC,2));
ContrastArray(1,MaskVoxels_spatial)=Contrast_FacePlaces_Image;

SimbasMask=reshape(SimbasArray,91,109,91);
ContrastImage=reshape(ContrastArray,91,109,91);
SimbasImage=make_nii(SimbasMask);
Contrast=make_nii(ContrastImage);
save_nii(Contrast,strcat(path2,'Contrast_FacesPlaces_NSB'));
%Clustering
C=SimbasImage.img;
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
NClusters125_NSB=length(IndexClusterSize_125);
NVoxels125_NSB=sum(IndexClusterSize_125);
NClusters64_NSB=length(IndexClusterSize_64);
NVoxels64_NSB=sum(IndexClusterSize_64);
NClusters27_NSB=length(IndexClusterSize_27);
NVoxels27_NSB=sum(IndexClusterSize_27);
BiggCl_NSB=IndexClusterSize_125(1);

save(strcat(path2,'Results_Table3_FacePlaces_task_NSB'),'IntWidth','NClusters125_NSB', 'NVoxels125_NSB','NClusters64_NSB','NVoxels64_NSB','NClusters27_NSB','NVoxels27_NSB','BiggCl_NSB')

%% Cluster image

for ii=1:NClusters64_CHSB
TesteImg=zeros(size(C));
Region=find(Aux1==IndexClusterSize_64(ii));
TesteImg(CC.PixelIdxList{Region(1)})=1;
BigCluster=make_nii(TesteImg);
save_nii(BigCluster,strcat(path1,'Results',sep,'Cluster',num2str(ii),'_NSB'))
end