% SHAPEL ------------------------------------------------------------------
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
function [SF,GDSF,DET] = SHAPEL(XI,ELXY)
% *************************************************************************
% COMPUTE SHAPE FUNCTIONS, DERIVATIVES AND DETERMINANT OF HEXAGONS
% *************************************************************************

XNODE = [...
    -1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
QUAR = 0.125;

SF = zeros(8,1);        DSF = zeros(3,8);

for I = 1:8
    XP = XNODE(1,I);
    YP = XNODE(2,I);
    ZP = XNODE(3,I);
    
    XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];
    
    SF(I) = QUAR*XI0(1)*XI0(2)*XI0(3);
    DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
    DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
    DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
end

GJ   = DSF*ELXY;
DET  = det(GJ);
GDSF = GJ\DSF;
end