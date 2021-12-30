
path_spm='/Users/michellemiranda/Documents/MATLAB/SPMTest/';

load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'HCP_Application',sep);
path2=strcat(path1,'Results',sep);


C3=load_nii(strcat(path_spm,'beta_0003.nii'));
CI3=C3.img;
C11=load_nii(strcat(path_spm,'beta_0011.nii'));
CI11=C11.img;

C7=load_nii(strcat(path_spm,'beta_0005.nii'));
CI7=C7.img;
C15=load_nii(strcat(path_spm,'beta_0013.nii'));
CI15=C15.img;

ContrastArray=(CI3+CI11)/2-(CI7+CI15)/2;
I2=load_nii(strcat(path,sep,'Talairach-labels-2mm.nii'));
C=I2.img;
[dimx,dimy,dimz]=size(C);
[iv,jv,kv]=ind2sub(size(C), find(C~=0) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);


Img_SPM_mask=zeros(size(C));
Img_SPM_mask(MaskVoxels)=ContrastArray(MaskVoxels);
save_nii(make_nii(Img_SPM_mask),strcat(path2,'Contrast_FacePlaces_SPM'));


SpmI=load_nii(strcat(path_spm,'spmT_0001.nii'));
Teste=SpmI.img;
TT=reshape(Teste,91*109*91,1);
MaskTT=find(abs(TT)>5.14);
TT2=zeros(size(TT));
TT2(MaskTT)=1;
C=reshape(TT2,91,109,91);

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
NClusters125_CHSB=length(IndexClusterSize_125);
NVoxels125_CHSB=sum(IndexClusterSize_125);
NClusters64_CHSB=length(IndexClusterSize_64);
NVoxels64_CHSB=sum(IndexClusterSize_64);
NClusters27_CHSB=length(IndexClusterSize_27);
NVoxels27_CHSB=sum(IndexClusterSize_27);
BiggCl_CHSB=ClusterSizes(1);

save(strcat(path2,'Results_Table3_FacePlaces_task_SPM'),'NClusters125_CHSB', 'NVoxels125_CHSB','NClusters64_CHSB','NVoxels64_CHSB','NClusters27_CHSB','NVoxels27_CHSB','BiggCl_CHSB')


