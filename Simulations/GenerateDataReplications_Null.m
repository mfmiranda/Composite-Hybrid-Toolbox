%This function will generate replications for the simulation studies
% 1) 100 task datasets with long term memory
% 2) 100 task datasets with short term memory
% 3) 100 null datasets with long term memory

load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'Simulations',sep);
path2=strcat(path1,'Replications_100',sep);
path3=strcat(path1,'Datasets',sep);

addpath(genpath(strcat(path,sep,'CodeDWT')))
addpath(genpath(strcat(path,sep,'RWFMM_Allcode')))
addpath(genpath(strcat(path,sep,'tools-for-nifti-and-analyse-image')))
load(strcat(path1,'X_Design'))
I2=load_nii(strcat(path,sep,'Talairach-labels-2mm.nii'));
I_RL=load_nii(strcat(path3,'Task_centered'));
C=I2.img(32:63,32:63,26:50);
C2=I_RL.img(32:63,32:63,26:50);
load(strcat(path1,'Sim_PhiMatrix'))
NN=100; %Total number of datasets

tic
for dataset=1:NN
    

K0=1;
K1=0.5;
K2=2.2;


 % SIGNAL GENERATION

[dimx,dimy,dimz]=size(C);


%Generate AR(2) process
NV=dimx*dimy*dimz;
Y=zeros(NV,NN+2);
for kk=1:NV
    for tt=1:NN
    Y(kk,tt+2)=0.9*Y(kk,tt+1)-0.8*Y(kk,tt)+normrnd(0,0.5);
    end
end

noise=Y(:,3:end);

% NN images with short-term noise
 % short-term noise
 for k=1:dimz    
    for i=1:dimx       
        for j=1:dimy  
            if k<dimz && i<dimx && j<dimy
               noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
               noise(i+1+(j-1)*dimx+(k-1)*dimx*dimy,:)+noise(i+j*dimx+(k-1)*dimx*dimy,:)+noise(i+(j-1)*dimx+k*dimx*dimy,:))/4;
            elseif k==dimz && i<dimx && j<dimy
               noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
               noise(i+1+(j-1)*dimx+(k-1)*dimx*dimy,:)+noise(i+j*dimx+(k-1)*dimx*dimy,:))/3;
            elseif k<dimz && i==dimx && j<dimy
               noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
               noise(i+j*dimx+(k-1)*dimx*dimy,:)+noise(i+(j-1)*dimx+k*dimx*dimy,:))/3;
            elseif k<dimz && i<dimx && j==dimy
               noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
               noise(i+1+(j-1)*dimx+(k-1)*dimx*dimy,:)+noise(i+(j-1)*dimx+k*dimx*dimy,:))/3;
           elseif k<dimz && i==dimx && j==dimy
             noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
                                              noise(i+(j-1)*dimx+k*dimx*dimy,:))/2;
           elseif k==dimz && i<dimx && j==dimy
             noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
                               noise(i+1+(j-1)*dimx+(k-1)*dimx*dimy,:))/2;
           elseif k==dimz && i==dimx && j<dimy  
             noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=(noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:)+...
                                              noise(i+j*dimx+(k-1)*dimx*dimy,:))/2;
            else
             noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,:)=noise(i+(j-1)*dimx+(k-1)*dimx*dimy,:);
            end           
      end        
    end    
 end

 %Long-term noise
 noise_long1=zeros(size(noise));
for rep=1:NN
    
for k=1:dimz    
    for i=1:dimx       
        for j=1:dimy 
         noise_long1(i+(j-1)*dimx+(k-1)*dimx*dimy,rep)=2*sin(2*pi*i/dimx)*noise(i+(j-1)*dimx+(k-1)*dimx*dimy,rep)+2*cos(2*pi*j/dimy)*noise(i+(j-1)*dimx+(k-1)*dimx*dimy,rep)+2*sin(2*pi*k/dimz)*noise(i+(j-1)*dimx+(k-1)*dimx*dimy,rep)+noise_short1(i+(j-1)*dimx+(k-1)*dimx*dimy,rep);
         end        
    end    
end

end 


for tt=1:NN
    ImgData(:,:,:,tt)=(double(C2))+K1*reshape(noise_long1(:,tt),dimx,dimy,dimz);
end 

% Running our proposed hybrid basis model

Y_Matrix=(reshape(ImgData,dimx*dimy*dimz,size(ImgData,4)))';
Y_MatrixC = bsxfun(@minus, Y_Matrix, mean(Y_Matrix));
filename=strcat(path2,'DataNull_LongTermNoise_',num2str(dataset));
save(filename,'Y_MatrixC')
end %end dataset generation