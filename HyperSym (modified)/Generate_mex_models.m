% Generate_mex_models -----------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function Generate_mex_models(MID,PROP,SRI)
% Remove warning from Matlab coder
warning('off','MATLAB:RMDIR:RemovedFromPath')
warning('off','MATLAB:DELETE:Permission')
% Make a copy of the '.m' file and rename it to 'fnc_mex.mexw64'
CopyGenerateRenameMove(MID,numel(PROP))
if SRI % For SRI implementation, create an isochoric and volumetric file
    CopyGenerateRenameMove([MID '_iso'],numel(PROP))
    CopyGenerateRenameMove([MID '_vol'],numel(PROP))
end
delete('fnc.m');
end
%% Operations with '.m' file
function CopyGenerateRenameMove(name,nPROP)
% Copy '.m' file in current folder, generate mex, and store it in folder 
% 'Functions'
copyfile(['Functions\' name '.m'],'fnc.m'); pause(.01);
Generate_mex(nPROP); pause(.01);
system(['rename ' 'fnc_mex.mexw64 ' name '.mexw64']); pause(.01);
movefile([name '.mexw64'],'Functions')
end
%% MEX file creation
function Generate_mex(nPROP)
% Configure configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
[cfg.GenerateReport,cfg.ReportPotentialDifferences]= deal(false);
cfg.InlineBetweenUserAndMathWorksFunctions = 'Speed';
% Create MEX file
disp('Compiling MEX files...');
addpath('Functions/');
% Define inputs of the MEX file: Right Cauchy-Green tensor and properties
C = coder.typeof(0,[3 3],[0 0]);
PROP = coder.typeof(0,[1 nPROP],[0 0]);
ARGS = {C; PROP}; % Store input arguments in cell array
cfg.Verbosity = 'Silent';
% Invoke MATLAB Coder
codegen -config cfg fnc -args ARGS % Invoke MATLAB Coder.
end