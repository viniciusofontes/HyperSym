% Strip -------------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
%
% References:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element analysis.
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

% Load increments [Start End Increment InitialFactor FinalFactor]
TIMS = [0.0 1.0 0.05 0.0 1.0]';

% Material ID
MID = 'MR';

% Material properties
psi2MPa = 6.8945757293E-3;

A10 = 24*psi2MPa;
A01 = 1.5*psi2MPa;
mu = 2*(A01 + A10);           % 2nd Lamé's parameter (MPa)
nu = 0.4;
K = 2*mu*(1+nu)/(3*(1-2*nu)); % Bulk modulus (MPa)
lambda = K - 2*mu/3;          % 1st Lamé's parameter (MPa)

PROP = [A10 A01 K];

% Set program parameters
ITRA = 30;    % maximum number of Newton-Raphson iterations'
NTOL = 6;     % maximum number of consecutive bisections
ATOL = 1E5;   % divergence value (if residual is greater than ATOL,
...             the bisection procedure is evoked)
TOL = 1E-6;   % tolerance convergence
    
N = [4 8 16 32 64]; % Mesh size             
[XYZ,LE,DIR] = deal(cell(length(N),1));  m = 1;

% SRI: Selects reduced integration (if available for the chosen model)
% false, full integration; true, selective reduced integration (SRI)
SRI = true;
    
while N(m) <= N(end)
    MeshSize = [N(m) N(m) 1];
    
    fprintf('Mesh #%d: %d x %d x %d\n',m,MeshSize(1),MeshSize(2),MeshSize(3))
    fprintf('----------------------------------------\n\n')
    
    [XYZ{m},LE{m}] = CreateMesh_Strip(BdBox,MeshSize);
    
    % Prescribed displacements [Node, DOF, Value]
    LeftNodes  = find(abs(XYZ{m}(:,1) - BdBox(1)) < eps);
    RightNodes = find(abs(XYZ{m}(:,1) - BdBox(2)) < eps); % non-zero displacement
    LowerNodes = find(abs(XYZ{m}(:,2) - BdBox(3)) < eps);
    BaseNodes  = find(abs(XYZ{m}(:,3) - BdBox(5)) < eps);
    
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
    
    % Contact elements
    LEC = [];
    
    MeshName = sprintf('_N=%d',MeshSize(1)*MeshSize(2)*MeshSize(3));
    
    % Mesh
    % PlotMesh_Strip({FileName; MeshName},XYZ{m},LE{m},SDISPT,EXTFORCE)
    
    DIR{m} = Directory({FileName; MeshName},MID,SRI);

    % Calling main function
    NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,MID,PROP,EXTFORCE,SDISPT,XYZ{m},LE{m},DIR{m},SRI)
    
    % Plotting Deformed Mesh
    ARGS = struct(...
    'axis',[0 Dim.L/2 + U 0 Dim.L/2 0 Dim.t],...
    'view',[55 25],...
    'xticks',0:30:ceil(Dim.L/2 + U),...
    'yticks',0:30:Dim.L/2,...
    'zticks',0:Dim.t);

    Plotter(XYZ{m},LE{m},DIR{m},ARGS)
    % StressField(DIR{m},{FileName; MeshName},MID,SRI,XYZ{m},LE{m},'Mises','CART')
    
    fprintf('\t\t *** Analysis converged! *** \t\n\n')

    if m == length(N)
        fprintf('\t *** Successful end of program! *** \t\n\n')
        break;
    else
        m = m + 1;
    end
end
%% Post-processing
PostProc_Strip(DIR,SRI,XYZ,LE,BdBox,N)