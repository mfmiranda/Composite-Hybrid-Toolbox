
%% Obtain data from HCP and save it into Y_MatrixC
I_RL=load_nii('tfMRI_WM_RL.nii');
I0_RL=I_RL.img;
[dimx,dimy,dimz,dimt]=size(I0_RL);
Y_Matrix=(double(reshape(I0_RL,dimx*dimy*dimz,dimt)))';
Y_MatrixC = bsxfun(@minus, Y_Matrix, mean(Y_Matrix));
save('ExampleSubjectData', 'Y_MatrixC', '-v7.3')
