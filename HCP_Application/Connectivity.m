load('ToolboxPath')
if Mac==1; sep='/'; else; sep='\'; end
path1=strcat(path,sep,'HCP_Application',sep);
path2=strcat(path1,'Results',sep);
addpath(genpath(strcat(path,sep,'tools-for-nifti-and-analyse-image')))
sq_RV=readmatrix(strcat(path2,'Rv2.txt'));
sq_RV_reordered=readmatrix(strcat(path2,'ConnectivityOrdered.txt'));

figure
imagesc(sq_RV_reordered)
title('Estimated Background Connectivity')
colormap('jet')
colorbar('FontSize',12)
saveas(gcf,strcat(path2,'ConnectivityOrdered.png'))

%Extract lower traingular matrix
LT=tril(sq_RV,-1);
test=reshape(LT,88804,1);
[testSort,ind]=sort(test,'descend');

k08=find(testSort>0.8); %corr>0.8
k07=find(testSort>0.7); %corr>0.8
%[testSort(k08) test(ind(k08))] %testing if they are equivalent
[row, col] = ind2sub(size(sq_RV), ind);
PairsSq=[row(k07), col(k07)]
%row and col give me the row and column in the RV matrix ordered by the
%connectivity
I2=load_nii(strcat(path,sep,'Talairach-labels-2mm.nii'));
C=I2.img;
[dimx, dimy,dimz]=size(C);
[valD numD] = howmany(C);
index_lesspart=find(numD>125);
Original_index=index_lesspart(2:end);
%Original_index gives me the index in the labels vectors
Labels=valD(Original_index);
TableLabels=[double(Labels(row(k07))),double(Labels(col(k07))),double(sq_RV(ind(k07)))];

%Plotting Selected Regions

SelectedRegions=unique([TableLabels(:,1);TableLabels(:,2)]);
NewMatrix_SqRV=zeros(size(SelectedRegions,1),size(SelectedRegions,1));

for i=1:size(SelectedRegions)
    for j=1:size(SelectedRegions)
     NewMatrix_SqRV(i,j)= sq_RV(find(Labels==SelectedRegions(i)),find(Labels==SelectedRegions(j)));
    end
end
tickpoints=[1:1:length(SelectedRegions)];
imagesc(NewMatrix_SqRV)
title('Estimated Background Connectivity - Selected Regions')
colormap('jet')
colorbar('FontSize',12)
yticklabels=num2str(SelectedRegions);
set(gca, 'XTick', tickpoints, 'XTickLabel', yticklabels) 
set(gca, 'YTick', tickpoints, 'YTickLabel', yticklabels) 
saveas(gcf,strcat(path2,'ConnectivityReordered_SelectedRegions.png'))

% NewMatrix_SqRV is the matrix that inputs in Brain View
% We choose a smaller subset of the 26 selected regions
% to facilitate visualization

NewSelect=[81,238,244,473,478,667,844,902,987];
%aux0_indexROI=find(SelectedRegions==NewSelect);
BrainView_SqRV=zeros(size(NewSelect,2),size(NewSelect,2));

for i=1:size(NewSelect,2)
    for j=1:size(NewSelect,2)
     BrainView_SqRV(i,j)= sq_RV(find(Labels==NewSelect(i)),find(Labels==NewSelect(j)));
    end
end

for k=1:size(NewSelect,2)
BrainView_SqRV(k,k)=-1;
end
writematrix(BrainView_SqRV, strcat(path2,"BrainViewRV.txt"),'Delimiter','tab');
