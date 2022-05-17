load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'Simulations',sep);

I2=load_nii(strcat(path,'Talairach-labels-2mm.nii'));
C=I2.img(32:63,32:63,26:50);
[valD numD] = howmany(C);
NVox=sum(numD);
NRoi=size(numD,1);

Phi_Matrix=zeros(NRoi,NVox);
 NN=100;
[dimx,dimy,dimz]=size(C);

for roi_i=1:NRoi

[iv,jv,kv]=ind2sub(size(C), find(C==valD(roi_i)) );
ssize=[dimx, dimy,dimz];
MaskVoxels=sub2ind(ssize,iv,jv,kv);
PhiArray=zeros(size(C));
PhiArray(MaskVoxels)=1;
Phi_Matrix(roi_i,:)=reshape(PhiArray,dimx*dimy*dimz,1);
end

index_lesspart=find(sum(Phi_Matrix,2)>125);%exclude small partitions
Phi_MatrixNew=Phi_Matrix(index_lesspart,:);
Phi_MatrixNew(1,:)=[];

 % SIGNAL GENERATION

Size_ROI=numD(index_lesspart(2:end));
AuxSize=cumsum(Size_ROI);
IndexPart=valD(index_lesspart);
IndexPart(1)=[];
save(strcat(path1,'Sim_PhiMatrix'),'Phi_MatrixNew','IndexPart','AuxSize')