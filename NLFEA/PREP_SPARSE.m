% PREP_SPARSE -------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
%
% Reference:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element Analysis.
%     Springer Science & Business Media, pp. 214-216.
% -------------------------------------------------------------------------
function [i,j] = PREP_SPARSE(LE)
% *************************************************************************
% PREPRE SPARSITY FOR THE STIFFNESS MATRIX
% *************************************************************************
[NE,NNE] = size(LE);              NDof = 3*NNE; % Number of Dofs per el.
[i,j] = deal(zeros(NE*NDof^2,1)); indexK = 0;

for el = 1:NE
    NKE = NDof^2;

    % Element DOFs
    eDof = reshape([3*LE(el,:)-2;3*LE(el,:)-1;3*LE(el,:)],NDof,1);
    
    % Triplets for stiffness matrix
    I = repmat(eDof ,1,NDof); J = I';
    i(indexK+1:indexK+NKE) = I(:);
    j(indexK+1:indexK+NKE) = J(:);
    indexK = indexK + NKE;
end
end