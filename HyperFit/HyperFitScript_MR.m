% HyperFitScript_MR -------------------------------------------------------
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
%
% Digitalized data from [1], with stress converted to MPa.

kg2MPa = 9.80665/10^2; % 1 kgf/cm² = 9.80665 N/(10 mm)²

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
%% Define the material model: Mooney-Rivlin
% neo-Hookean
syms A10 A01 lambda1 lambda2 lambda3
I1 = lambda1^2 + lambda2^2 + lambda3^2;
I2 = lambda1^2*lambda2^2 + lambda2^2*lambda3^2 + lambda3^2*lambda1^2;
W = A10*(I1 - 3) + A01*(I2 - 3);
PROP = [A10 A01]; % Array of properties

% Use only a portion of the data for the data fit
portion = input('Wanna fit to the whole data? ''Y''/''N'' [Y]\n','s');
fprintf('\n')

if isempty(portion), prompt = 'Y'; end

if ismember(portion,{'N','n'})
    splice = @(x,val) x(x(:,1) <= val,:); xval = 2.0;

    for i = 1:length(Test)
        Test(i).data = splice(Test(i).data,xval);
    end
end
%% Call HyperFit
opt = struct(...
    'N',1,...     Number of optimization runs
    'pL',[],...   Lower boundary of p
    'pU',[],...   Upper boundary of p
    'p0',[],...   First initial guess
    'cc',true ... Check for consistency conditions
    );
% HyperFit returns fitted parameters (p), and error metrics (Error, res)
[p,Error,res] = HyperFit(W,PROP,Test,opt); 
%% Plot results
% Plot fitted model's results alongside neoHookean's
HyperFit_Comparison_neoHookean

% Plot fitted model's relative error
HyperFit_Relative_Error

if ismember(portion,{'Y','y'})
    % Adjust plot with all tests for all tests: YLim up to 8MPa
    ax = figure(1).Children(end); ax.YLim = [0 8];
    % Adjust plot with all tests for all tests (rel. err.): YLim up to 100%
    ax = figure(2).Children(end); ax.YLim = [0 100];
else
    % Adjust plot with all tests to display only for 1 <= lambda <= 2
    for i = 8:-2:2
        ax = figure(1).Children(i); ax.XLim = [1 2];
    end
    ax = figure(1).Children(8); ax.YLim = [0 1.2];

    % Extend the dashed line in all tests (rel. error) plot to lambda = 1
    ax = figure(2).Children(end); ax.XLim = [1 2]; ax.YLim = [0 20];
    [ax.Children(1).XData,ax.Children(1).YData] = deal([1 2],[5 5]);
end