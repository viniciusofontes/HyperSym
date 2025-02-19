% Uniaxial_Traction_Comp_Models -------------------------------------------
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
%% Pre-processing
clc, clear, close all

% Go to the directory of this file
FileName = mfilename;
cd(fileparts(which(FileName)))

addpath(genpath('./')); % Use all folders and subfolders

set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')

set(0,'DefaultFigureWindowStyle','normal')
%% Problem definition
% Nodal coordinates
XYZ = [...
    0 -1 0;
    1 -1 0;
    0  0 0;
    1  0 0;
    0 -1 1;
    1 -1 1;
    0  0 1;
    1  0 1];

% Element connectivity
LE = [1 2 4 3 5 6 8 7];

% External forces [Node DOF Value]
EXTFORCE = [];

% Prescribed displacements [Node DOF Value]
SDISPT = [...
    1 1 0;1 2 0;1 3 0;
    2 1 1;2 2 0;2 3 0;
    3 1 0;3 2 0;3 3 0;
    4 1 1;4 2 0;4 3 0;
    5 1 0;5 2 0;5 3 0;
    6 1 1;6 2 0;6 3 0;
    7 1 0;7 2 0;7 3 0;
    8 1 1;8 2 0;8 3 0];

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS = [0 1.0 0.05 0.0 1.0]';

% Contact elements
LEC = [];

% Material ID
MID = {'SVK'; 'mSVK1'; 'mSVK2'; 'mSVK3'; 'nH1'; 'nH2'; 'nH3'};
nMAT = length(MID);            DIR = cell(nMAT,1);

% Material properties
E = 1000;
nu = 0.3;
mu = E/(2*(1 + nu));
lambda = E*nu/((1 + nu)*(1 - 2*nu));
K = lambda + 2/3*mu;

% Set program parameters
ITRA = 30;    % maximum number of Newton-Raphson iterations'
NTOL = 6;     % maximum number of consecutive bisections
ATOL = 1E5;   % divergence value (if residual is greater than ATOL,
...             the bisection procedure is evoked)
TOL = 1E-6;   % tolerance convergence

% SRI: Selects reduced integration (if available for the chosen model)
% false: full integration; true: selective reduced integration (SRI)
SRI = false;
% nMAT=1;
for m = 1:nMAT
    fprintf('Model #%g: %s\n',m,MID{m})
    fprintf('----------------------------------------\n\n')
       
    if ismember(MID{m},{'SVK'; 'mSVK1'; 'mSVK2'; 'mSVK3'; 'nH1'; 'nH2'; 'nH3'}) % SVK or nH
        PROP = [lambda mu];
    elseif strcmp(MID{m},'MR') % Mooney-Rivlin
        PROP = [A10 A01 K];
    else
        errordlg('Wrong material ID!','Input Error'); return;
    end
    
    DIR{m} = Directory({FileName},MID{m},SRI);
    
    % Calling main function
    NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,MID{m},PROP,EXTFORCE,SDISPT,XYZ,LE,DIR{m},SRI);

    % Plotting Deformed Mesh
    ARGS = struct(...
        'axis',[0 2 -1 0 0 1],...
        'view',[45 25],...
        'xticks',0:0.2:2);
    
    if m == nMAT, Plotter(XYZ,LE,DIR{m},ARGS); end
    
    fprintf('  *** Material analysis has converged! ***  \n\n')
end

fprintf('\n\t ***  Successful end of program! *** \t\n\n')
%% Post-processing
PostProc_UniaxialTraction_Comp_Models(DIR,SRI,XYZ,LE,8,MID,E,nu)