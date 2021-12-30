function W=GetW_Michelle(model,D)

%%%%% W=GetW(model,D)
%%%%%   Compute cross-products matrices
%%%%%           X'X       X'Z       X'D
%%%%%           Z'X       Z'Z       Z'D
%%%%%           D'X       D'Z       D'D
% Input:    
%
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%           D = (n x K) matrix of wavelet coefficients
%
% Output:   W = structure with elements:
%           XtX,XtZ,XtD,ZtD,DtD

W.XtX=model.X'*model.X;
W.XtD=model.X'*D;
W.DtD=D'*D;
