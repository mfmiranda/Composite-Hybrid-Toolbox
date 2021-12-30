% This script will prepare the data and run the transform
% it will then split the data into 10 different files to be run
% independently
load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'HCP_Application',sep);
path2=strcat(path1,'Results',sep);
%path3=strcat(path1,'Datasets',sep);

addpath(genpath(strcat(path,sep,'CodeDWT')))
addpath(genpath(strcat(path,sep,'RWFMM_Allcode')))
addpath(genpath(strcat(path,sep,'tools-for-nifti-and-analyse-image')))

load(strcat(path1,'DesignMatrix_100307.mat'))
load(strcat(path1,'ExampleSubjectData.mat'))

% No Spatial Basis - we will look at the same voxels we started for the
% other methods


Auxiliar=sum(Y_MatrixC,1);
MaskVoxels_spatial=find(Auxiliar~=0);
Y_MatrixNonZeros=Y_MatrixC(:,MaskVoxels_spatial);
    
%% Wavelet Transform on time domain
Y_star=Y_MatrixNonZeros';
wavespecs0.wavelet='db1';
wavespecs0.ndim=1;
wavespecs0.compress=1;
padding=107;
Y0=[Y_star,repmat(0,size(Y_star,1),padding)]; 
%Add padding until numer of basis is even. This step makes sure that Y=XB+E is the same as
%WY=WXB+WE

[Dtrans,wavespecs]=DWT_rows(Y0,wavespecs0);

XCov = bsxfun(@minus, DesignMatrix, mean(DesignMatrix));
XCov2 = bsxfun(@rdivide, XCov, std(DesignMatrix));
Xt=[XCov2;repmat(0,padding,size(XCov2,2))];

[DX,wavespecsX]=DWT_rows(Xt',wavespecs);
save('DX','DX')

kk=size(Y_MatrixNonZeros,2)/10;
for ii=1:10
D1=Dtrans((ii-1)*kk+1:ii*kk,:);
save(strcat(path1,'DTransform',num2str(ii)),'D1','wavespecs')
end


