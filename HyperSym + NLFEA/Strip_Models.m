% Strip_Models ------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% 
% References:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element Analysis.
%     Springer Science & Business Media, pp. 214-216.
% [2] Oden, J.T. (1972). Finite Elements of Nonlinear Continua, Courier 
%     Corporation, pp. 304-308.
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
    'L',8*pol2mm,...    % length (mm)
    't',0.05*pol2mm ... % thickness (mm)
    );

BdBox = [0 Dim.L/2 0 Dim.L/2 0 Dim.t];

% Prescribed displacement (mm) (Ogden, 1972)
stretch = 3; U = Dim.L*(stretch - 1)/2;

% Nodal coordinates and element connectivity
MeshSize = [16 16 1];
% MeshSize = [2 2 1];
[XYZ,LE] = CreateMesh_Strip(BdBox,MeshSize);

% Prescribed displacements [Node, DOF, Value]
LeftNodes  = find(abs(XYZ(:,1) - BdBox(1)) < eps);
RightNodes = find(abs(XYZ(:,1) - BdBox(2)) < eps); % non-zero displacement
LowerNodes = find(abs(XYZ(:,2) - BdBox(3)) < eps);
BaseNodes  = find(abs(XYZ(:,3) - BdBox(5)) < eps);

SDISPT = [...
    LeftNodes     ones(length(LeftNodes),1)   zeros(length(LeftNodes),1);
    RightNodes   ones(length(RightNodes),1) U*ones(length(RightNodes),1);
    RightNodes 2*ones(length(RightNodes),1)  zeros(length(RightNodes),1);
    LowerNodes 2*ones(length(LowerNodes),1)  zeros(length(LowerNodes),1);
    BaseNodes   3*ones(length(BaseNodes),1)   zeros(length(BaseNodes),1)];

% Sorting and eliminating duplicates dofs according to node number
SDISPT = sortrows(unique(SDISPT,'rows'));

% External forces [Node, DOF, Value]
EXTFORCE = [];

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS = [0.0 1.0 0.05 0.0 1.0]';

% Contact elements
LEC = [];

% Mesh
PlotMesh_Strip({FileName},XYZ,LE,SDISPT,EXTFORCE)

% Material ID
MID = {'SVK';'mSVK1';'mSVK2';'mSVK3';'nH1';'nH2';'nH3';'MR'};
nMAT = length(MID);            DIR = cell(nMAT,1);

% Material properties
psi2MPa = 6.8945757293E-3;

A10 = 24*psi2MPa;
A01 = 1.5*psi2MPa;
mu = 2*(A01 + A10);           % 2nd Lamé's parameter (MPa)
nu = 0.4;
K = 2*mu*(1+nu)/(3*(1-2*nu)); % Bulk modulus (MPa)
lambda = K - 2*mu/3;          % 1st Lamé's parameter (MPa)

% Set program parameters
ITRA = 30;    % maximum number of Newton-Raphson iterations'
NTOL = 6;     % maximum number of consecutive bisections
ATOL = 1E5;   % divergence value (if residual is greater than ATOL,
...             the bisection procedure is evoked)
TOL = 1E-6;   % tolerance convergence

% SRI: Selects reduced integration (if available for the chosen model)
% false, full integration; true, selective reduced integration (SRI)
SRI = false;
runFEA = true; % true: Run FEA; false: Only view results 
plotMesh = true; % true: Plot deformed mesh (animation); false: Don't
for m = 1:nMAT
    fprintf('Model %g: %s\n',m,MID{m})
    fprintf('----------------------------------------\n\n')
       
    if ismember(MID{m},{'SVK'; 'mSVK1'; 'mSVK2'; 'mSVK3'; 'nH1'; 'nH2'; 'nH3'}) % SVK or nH
        PROP = [lambda mu];
    elseif strcmp(MID{m},'MR') % Mooney-Rivlin
        PROP = [A10 A01 K];
        SRI = true;
    else
        errordlg('Wrong material ID!','Input Error'); return;
    end
    
    DIR{m} = Directory({FileName},MID{m},SRI);
    
    % Calling main function    
    if runFEA
        NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,MID{m},PROP,EXTFORCE,SDISPT,XYZ,LE,DIR{m},SRI)
    end
    
    if plotMesh
        % Plotting Deformed Mesh
        ARGS = struct(...
            'axis',[0 Dim.L/2 + U 0 Dim.L/2 0 Dim.t],...
            'view',[55 25],...
            'xticks',0:30:ceil(Dim.L/2 + U),...
            'yticks',0:30:Dim.L/2,...
            'zticks',0:Dim.t);
        Plotter(XYZ,LE,DIR{m},ARGS)
    end

    fprintf('\t\t *** Analysis converged! *** \t\n\n')
end
%% Post-processing
PostProc_Strip_Models(DIR,SRI,XYZ,LE,BdBox,MID,'2D')