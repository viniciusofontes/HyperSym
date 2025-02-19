% Generate_mex_models -----------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
%% MEX file creation function
clc, clear, close all
% Go to this file's folder
cd(fileparts(which(mfilename)));

% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.InlineBetweenUserAndMathWorksFunctions = "Speed";
% cfg.EnableJIT = true;
% Create MEX file
disp('Compiling MEX files...');

% MID
MID = coder.newtype('char', [1 20], [0 1]);
% PROP
PROP = coder.typeof(0,[1 Inf],[0 1]);
% UPDATE,LTAN
[UPDATE,LTAN] = deal(coder.typeof(true));
% XYZ
XYZ = coder.typeof(0,[Inf 3],[1 0]);
% LE
LE = coder.typeof(0,[Inf 8],[1 0]);
% DISPTD,FORCE,INTFORCE
DISPTD = coder.typeof(0,[Inf 1],[1 0]);
% DSFS
DSFS = coder.typeof(0,[3 8 Inf],[0 0 1]);

ARGS = {...
    MID;
    PROP;
    UPDATE;
    LTAN;
    XYZ;
    LE;
    DISPTD;
    DSFS...    
};

% Invoke MATLAB Coder.
codegen -config cfg HYPER3D -args ARGS
codegen -config cfg HYPER3D_SRI -args ARGS
% Rename the MEX files to match their non-MEX counterparts
system(['rename ' 'HYPER3D_mex.mexw64 HYPER3D.mexw64']);
system(['rename ' 'HYPER3D_SRI_mex.mexw64 HYPER3D_SRI.mexw64']);
