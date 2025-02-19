% HyperSym ----------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function [Sv,Dv,C] = HyperSym(W)
%% Symbolic Variables
syms I1 I2 I3 J J1 J2 real
C = sym('C', [3,3],'real'); % Right Cauchy–Green def. tensor, C = F'*F
%% Principal Invariants of C
mI1 = trace(C);              mI2 = (trace(C)^2 - trace(C^2))/2;
mJ = sqrt(det(C));           mI3 = mJ^2;
%% Reduced Invariants of C
mJ1 = mI3^(-1/3)*mI1;        mJ2 = mI3^(-2/3)*mI2;
%% Strain Energy Function, W(C)
W = subs(W,{I1,I2,I3,J,J1,J2},{mI1,mI2,mI3,mJ,mJ1,mJ2});
%% Derivation of PK2 Stress Tensor, S
S = sym(zeros(3,3));
for i = 1:3
    for j = 1:3
        S(i,j) = 2*diff(W(1),C(i,j));
    end
end
S = (S + transpose(S))/2; % Enforce symmetry of S
%% Derivation of Elasticity Tensor, D
D = sym(zeros(3,3,3,3));
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                D(i,j,k,l) = 2*diff(S(i,j),C(k,l));
            end
        end
    end
end
%% Impose symmetry and write in Voigt notation
S = simplify(subs(S,[C(2,1),C(3,2),C(3,1)],[C(1,2),C(2,3),C(1,3)]));
D = simplify(subs(D,[C(2,1),C(3,2),C(3,1)],[C(1,2),C(2,3),C(1,3)]));
Sv = [S(1,1) S(2,2) S(3,3) S(1,2) S(2,3) S(1,3)].';
Dv = [...
    D(1,1,1,1) D(1,1,2,2) D(1,1,3,3) D(1,1,1,2) D(1,1,2,3) D(1,1,3,1);
    D(2,2,1,1) D(2,2,2,2) D(2,2,3,3) D(2,2,1,2) D(2,2,2,3) D(2,2,3,1);
    D(3,3,1,1) D(3,3,2,2) D(3,3,3,3) D(3,3,1,2) D(3,3,2,3) D(3,3,3,1);
    D(1,2,1,1) D(1,2,2,2) D(1,2,3,3) D(1,2,1,2) D(1,2,2,3) D(1,2,3,1);
    D(2,3,1,1) D(2,3,2,2) D(2,3,3,3) D(2,3,1,2) D(2,3,2,3) D(2,3,3,1);
    D(3,1,1,1) D(3,1,2,2) D(3,1,3,3) D(3,1,1,2) D(3,1,2,3) D(3,1,3,1)];
end