% Thick_Walled_Cylinder_Incomp_Models -------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% 
% References:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element Analysis.
%     Springer Science & Business Media, pp. 214-216.
% [2] ANSYS, Inc. (2021). Ansys Mechanical APDL Verification Manual, 
%     pp. 157-158. Example VM56: Hyperelastic Thick Cylinder Under Internal
%     Pressure.
% [3] Oden, J.T. (1972). Finite Elements of Nonlinear Continua, Courier 
%     Corporation, pp. 325-331.
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
% Dimensions
pol2mm = 25.4;

Dim = struct(...        
    'ri',7*pol2mm,...        % inner radious (mm)
    'ro',18.625*pol2mm,...   % outer radious (mm)
    'h',0.775*pol2mm,...     % section height (mm)
    'alpha',90 ...           % section angle (degrees)
    );

BdBox = [Dim.ri Dim.ro Dim.ri Dim.ro 0 Dim.h];

% Nodal coordinates and element connectivity
MeshInfo = [20 20 1];

% Loads
psi2MPa = 6.8945757293E-3;
PRESS.Pi = 150*psi2MPa;      % inner pressure (MPa)
PRESS.EREF = MeshInfo(1);

[XYZ,LE] = CreateMesh_Cylinder(Dim,BdBox,MeshInfo);

% Prescribed displacements [Node, DOF, Value]
LeftNodes  = find(abs(XYZ(:,1)) < eps);
LowerNodes = find(abs(XYZ(:,2)) < eps);
AllNodes   = (1:length(XYZ))';

SDISPT = [...
     LeftNodes    ones(length(LeftNodes),1)  zeros(length(LeftNodes),1);
    LowerNodes 2*ones(length(LowerNodes),1) zeros(length(LowerNodes),1);
      AllNodes   3*ones(length(AllNodes),1)   zeros(length(AllNodes),1)];

% Sorting and eliminating duplicates dofs according to node number
SDISPT = sortrows(unique(SDISPT,'rows'));
    
% External forces (pressure magnitude)
EXTFORCE = PressureLoad(XYZ,LE,MeshInfo(1),PRESS.Pi);
EXTFORCE = sortrows(EXTFORCE);

% Mesh
PlotMesh_Cylinder(FileName,XYZ,LE,SDISPT,EXTFORCE);

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS = [0.0 1.0 0.1 0.0 1.0]';

% Contact elements
LEC = [];

% Material ID
% MID = {'MR'; 'Yeoh'};
MID = {'MR'};
nMAT = length(MID);    DIR = cell(nMAT,1);

% Material properties for MR
A10 = 80*psi2MPa;
A01 = 20*psi2MPa;
E = 6*(A10 + A01);
nu = 0.49;
K = E/(3*(1 - 2*nu));
mu = 2*(A10 + A01);    % See Kim p. 189
lambda = K - 2*mu/3;

A10Y = 0.645; A20Y = 0.0491; A30Y = -4.15E-4;
D1 = 2/K;     D2 = 1E9;      D3 = D2;

% Set program parameters
ITRA = 30;    % maximum number of Newton-Raphson iterations'
NTOL = 6;     % maximum number of consecutive bisections
ATOL = 1E5;   % divergence value (if residual is greater than ATOL,
...             the bisection procedure is evoked)
TOL = 1E-2;   % tolerance convergence

% SRI: Selects reduced integration (if available for the chosen model)
% false, full integration; true, selective reduced integration (SRI)
SRI = false;

for m = 1:nMAT
    fprintf('Model #%g: %s\n',m,MID{m})
    fprintf('----------------------------------------\n\n')
    
    if strcmp(MID{m},'MR') % Mooney-Rivlin
        PROP = [A10 A01 K];
    elseif strcmp(MID{m},'Yeoh')
        PROP = [A10Y A20Y A30Y D1 D2 D3];
    else
        errordlg('Wrong material ID!','Input Error'); return
    end
    
    DIR{m} = Directory({FileName},MID{m},SRI);
    
    % Calling main function
    NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,MID{m},PROP,EXTFORCE,SDISPT,XYZ,LE,DIR{m},SRI,PRESS);
    
    % Plotting Deformed Mesh
    ARGS = struct(...
        'axis',[0 600 0 600 0 ceil(Dim.h)],...
        'view',[40 30],...
        'zticks',[0 ceil(Dim.h)]);
    
    Plotter(XYZ,LE,DIR{m},ARGS)
    
    fprintf('\n  *** Material analysis has converged! ***  \n\n')
end
%% Post-processing
PostProc_Cylinder_Incomp_Models(DIR,SRI,XYZ,LE,BdBox,MID)