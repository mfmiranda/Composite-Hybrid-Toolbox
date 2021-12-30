%filename_data is the path of the folder where the downloaded data is
filename_data='/Users/michellemiranda/Documents/100307/MNINonLinear/Results/tfMRI_WM_RL/';
load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'HCP_Application',sep);

T=405; %number of time frames
p=8; %number of coefficients

%Reading EV's & creating the boxcar function

stringEV=cell(p,1);
stringEV{1}='2bk_body.txt';
stringEV{2}='2bk_faces.txt';
stringEV{3}='2bk_places.txt';
stringEV{4}='2bk_tools.txt';
stringEV{5}='0bk_body.txt';
stringEV{6}='0bk_faces.txt';
stringEV{7}='0bk_places.txt';
stringEV{8}='0bk_tools.txt';

% HRF
params.x = 6;
params.y = 16;
params.z = 1/6;
HRF = @(t) gampdf(t, params.x, 1) - params.z*gampdf(t, params.y, 1);

BoxStimulus=zeros(T,p); %Boxcar functions
dt_aux=0.36:0.72:301; %Time where frames were recorded

% Compute derivative
HRF_der_dt_aux=zeros(T,1);
HRF_der = @(x) (x^4*exp(-x))/24 - (x^5*exp(-x))/120 - (x^14*exp(-x))/523069747200 + (x^15*exp(-x))/7846046208000;
for i=1:T
    HRF_der_dt_aux(i,1)=HRF_der(dt_aux(i));
end

for i=1:p
filenameEV=stringEV{i};
EVFile=fopen(strcat(filename_data,'EVs',sep,filenameEV),'r');
A=fscanf(EVFile,'%f',3);
[index_frames]= FunctionEVs(A(1),A(2),0.72);
BoxStimulus(index_frames,i)=1;
end
%Convolution

B=zeros(T,p);
B2=zeros(T,p);
for i=1:8
    aux_hrf=conv(HRF(dt_aux),BoxStimulus(:,i));
    aux_dhrf=conv(HRF_der_dt_aux,BoxStimulus(:,i));
B(:,i)=aux_hrf(1:405);
B2(:,i)=aux_dhrf(1:405);
end

DesignMatrix=zeros(T,2*p);
for j=1:8
    DesignMatrix(:,(j-1)*2+1)=B(:,j);
    DesignMatrix(:,(j-1)*2+2)=B2(:,j);
end

save(strcat(path1,'DesignMatrix_100307'),'DesignMatrix')
imagesc(DesignMatrix)
colormap('jet')

