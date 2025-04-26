% ModHyperSymScript -------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
%
% Description: This script defines the hyperelastic model properties and
% strain energy density function (W) needed to derive the expressions for
% the 2nd Piola-Kirchhoff stress tensor (S) and the tangent elasticity
% tensor (D) by calling "HyperSym.m". These tensors are computed from the
% components of the right Cauchy-Green strain tensor (C). The quantities
% are converted to Voigt notation (Cv,Sv,Dv), which is the convention
% adopted in the generated functions.
%
% References:
% [1] Curnier, A. (1994). Computational Methods in Solid Mechanics.
%     Springer Science & Business Media
% [2] Klarbring, A. & Strömberg, N. (2013). Topology optimization of
%     hyperelastic bodies including non-zero prescribed displacements.
%     Structural and Multidisciplinary Optimization, 47(1), 37-48
% -------------------------------------------------------------------------
%% Pre-processing
clc; clear; close all
generateMEX = true; % Generate mex version of the  functions?
cd(fileparts(which(mfilename))); % Go to this file's folder
t = tic;
%% Symbolic variables
syms I1 I2 J J1 J2 lambda mu A10 A20 A30 A01 K D1 D2 D3 real
% Cell array to hold the strain energy density functions
MID = {'SVK'; 'mSVK1'; 'mSVK2'; 'mSVK3'; 'nH1'; 'nH2'; 'nH3'; 'MR'; 'Yeoh'};
W = cell(length(MID),1);
PROP = cell(length(MID),1);       % Array of properties for each model
[PROP{1:7}] = deal([lambda mu]);  % SVK and nH-based models
PROP{8} = [A10 A01 K];            % MR
PROP{9} = [A10 A20 A30 D1 D2 D3]; % Yeoh 3 parameters
%% Compressible models, i.e. W = W(I1,I2,J)
% Materials 1-4: St. Venant-Kirchhoff
W{1} = 1/8*lambda*(I1 - 3)^2 + mu/4*(I1^2 - 2*I2 - 2*I1 + 3);
W{2} = lambda/2*(log(J))^2  + mu/4*(I1^2 - 2*I2 - 2*I1 + 3);
W{3} = lambda*(J - log(J) - 1) + mu/4*(I1^2 - 2*I2 - 2*I1 + 3);
W{4} = lambda/2*(J - 1)^2 + mu/4*(I1^2 - 2*I2 - 2*I1 + 3);
% Materials 5-7: modified Neo-Hookean
W{5} = lambda/2*(log(J))^2 + mu/2*(I1 - 3) - mu*log(J);
W{6} = lambda*(J - log(J) - 1) + mu/2*(I1 - 3) - mu*log(J);
W{7} = lambda/2*(J - 1)^2 + mu/2*(I1 - 3) - mu*log(J);
%% Nearly incompressible models, i.e. W = Wiso(J1,J2) + Wvol(J)
% Wiso: iscochoric term, should be stored as the 1st component of W
% Wvol: volumetric term, should be stored as the 2nd component of W
W{8} = [...
    A10*(J1 - 3) + A01*(J2 - 3),... % Wiso
    K/2*(J - 1)^2];                 % Wvol
W{9} = [...
    A10*(J1 - 3) + A20*(J1 - 3)^2 + A30*(J1 - 3)^2, ... % Wiso
    1/D1*(J - 1)^2 + 1/D2*(J - 1)^4 + 1/D3*(J - 1)^6];  % Wvol
%% Genrate and store hyperelastic function
fprintf('Generating models...\n\n')
Dir = [cd '/Functions/'];  % Directory to save the functions
Controls = struct('File',[Dir 'fnc'],'Optimize',true);
if ~exist(Dir,'dir'), mkdir(Dir); end
for m = 1:length(MID)
    tic    
    fprintf('Model #%g: %s\n',m,MID{m})
    % Evaluate stress and Elasticity tensors        
    if length(W{m}) == 1 % Models without isochoric/volumetric split 
        [Sv,Dv,C] = HyperSym(W{m});        
        % Store tensors as column vector (for MEX compatibility)
        Dv = Dv(:); 
        % Generate Matlab functions
        Controls.File = [Dir MID{m}];
        matlabFunction(Sv,Dv,'Vars',{C,PROP{m}},Controls);
        % Generate MEX function
        if generateMEX,Generate_mex_models(MID{m},PROP{m},false); end
    else % Models with isochoric/volumetric split   
        [Siso,Diso,C] = HyperSym(W{m}(1));
        [Svol,Dvol] = HyperSym(W{m}(2));
        % Store tensors as column vector (for MEX compatibility)
        Sv = simplify(Siso+Svol);
        Dv = simplify(Diso+Dvol);
        Dv = Dv(:); Diso = Diso(:); Dvol = Dvol(:);
        % Generate Matlab functions
        Controls.File = [Dir MID{m}];
        matlabFunction(Sv,Dv,'Vars',{C,PROP{m}},Controls);
        % Generate Matlab functions (only isochoric part)
        Controls.File = [Dir MID{m} '_iso'];
        matlabFunction(Siso,Diso,'Vars',{C,PROP{m}},Controls);
        % Generate Matlab functions (only volumetric part)
        Controls.File = [Dir MID{m} '_vol'];
        matlabFunction(Svol,Dvol,'Vars',{C,PROP{m}},Controls);
        % Generate MEX function
        if generateMEX,Generate_mex_models(MID{m},PROP{m},true); end
    end
    fprintf('Time elapsed: %1.1f seconds.\n\n',toc)
end
%% Output information
fprintf('Done!\n\n')
if     toc(t) <= 60, fprintf('Total time: %1.1f seconds. \n',toc(t))
elseif toc(t) <= 3600, fprintf('Total time: %1.1f minutes. \n',toc(t)/60)
end