% EvalStresses ------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
%% Pre-processing
clc; clear; close all
%% Symbolic variables
syms lambda mu gamma
assume(gamma,'real')
%% Compute kinematic variables
I = sym(eye(3));        % Identity tensor
F = I; F(1,1) = gamma;  % Deformation gradient
C = F'*F;               % Right Cauchy-Green strain
Cinv = inv(C);          % Inverse of C
E = (C - I)/2;          % Green-Lagrange strain tensor
J = det(F);             % Jacobian of deformation
%% Analytical expressions for PK2 (Klarbring & Strömberg, 2013)
PK2{1} = [lambda*trace(E)*I 2*mu*E];
PK2{2} = [lambda*log(J)./C 2*mu*E];
PK2{3} = [lambda*(J - 1)./C 2*mu*E];
PK2{4} = [lambda*(J - 1)*J./C 2*mu*E];
PK2{5} = [lambda*log(J)./C  mu*(I - Cinv)];
PK2{6} = [lambda*(J - 1)./C  mu*(I - Cinv)];
PK2{7} = [lambda*(J - 1)*J./C  mu*(I - Cinv)];
%% Iterate through models and push-forward PK2
sigma11 = sym(zeros(7,2));

for m = 1:7
    cauchyL = 1/J*F*PK2{m}(:,1:3)*F';  cauchyR = 1/J*F*PK2{m}(:,4:6)*F';
    sigma11(m,:) = simplify([cauchyL(1,1) cauchyR(1,1)]);
    
    fprintf('sigma_%d(1,1) = \n',m); pretty((sigma11(m,:)))
end