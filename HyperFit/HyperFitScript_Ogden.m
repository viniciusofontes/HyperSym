% HyperFitScript_Ogden ----------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
%
% Reference:
% [1] Treloar, L.R.G. (1944). Stress-strain data for vulcanized rubber 
%     under various types of deformation. Rubber Chemistry and Technology, 
%     17(4), 813-825.
% -------------------------------------------------------------------------
%% Pre-processing
clc; clear; close all
addpath(genpath('./')) % Use all folders and subfolders
rng('default')         % Initial guesses
%% Input data from experimental tests (nominal stress vs principal stretch)
% Test: Cell array with structs containing 'name' and 'data' of the tests
% Test{i}.data: matrix containing principal stretch (1st collumn) and
% nominal stress (2nd collumn). 
% All test data must have same unit!
kg2MPa = 9.80665/10^2; % 1 kgf/cm^2 = 9.80665 N/(10 mm)^2
% Digitalized data from [1], with stress converted to MPa.
% ua: Simple Elongation (uniaxial stress)
ua_exp = readmatrix('Treloar Data/Treloar_uniaxial.csv');
ua_exp(:,2) = ua_exp(:,2)*kg2MPa; % Converting units
Test(1) = struct(...
    'name','uniaxial',...
    'data',ua_exp);
% eb: 2-Dimensional Extension (equi-biaxial stress)
eb_exp = readmatrix('Treloar Data/Treloar_equibiaxial.csv');
eb_exp(:,2) = eb_exp(:,1).*eb_exp(:,2)*kg2MPa; % Converting units
Test(2) = struct(...
    'name','equibiaxial',...
    'data',eb_exp(:,1:2));
% ps: Pure Shear
ps_exp = readmatrix('Treloar Data/Treloar_pureshear.csv');
ps_exp(:,2) = ps_exp(:,2)*kg2MPa;
Test(3) = struct(...
    'name','pureshear',...
    'data',ps_exp);
%% Define the material model: Ogden
syms lambda1 lambda2 lambda3
M = 4; % Number of pairs alpha and mu
W = sym(0);                           
alpha = sym('alpha',[M 1],'real');    mu = sym('mu',[M 1],'real');
for i = 1:M
    W = W + (mu(i)/alpha(i))*(...
        lambda1^alpha(i) + lambda2^alpha(i) + lambda3^alpha(i) - 3);
end
% Array of properties
PROP = reshape([alpha,mu]',1,2*M); % PROP = [alpha1, mu1, alpha2, mu2, ...]
%% Call HyperFit
opt = struct(...
    'N',5,...     % Number of optimization runs
    'pL',-10,...  % Lower boundary of p
    'pU',+10,...  % Upper boundary of p
    'p0',[],...   % First initial guess
    'cc',true ... % Check for consistency conditions
    );
% HyperFit returns fitted parameters (p), and error metrics (Error, res)
[p,Error,res] = HyperFit(W,PROP,Test,opt); 
%% Plot results
% Plot fitted model's results alongside neoHookean's
HyperFit_Comparison_neoHookean
% Plot fitted model's relative error
HyperFit_Relative_Error 
ax = figure(2).Children(end); ax.YLim = [0 10];