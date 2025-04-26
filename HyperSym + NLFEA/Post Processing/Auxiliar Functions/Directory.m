% Directory ---------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function DIR = Directory(FileName,MID,SRI)
% Create specif folders.

Folder = replace(FileName{1},'_',' '); % Save the output on a folder without '_' 
Path = [pwd '/Output/' Folder];        % Define the output folder

if ~exist([Path '/.mat Files'],'dir'), mkdir([Path '/.mat Files']); end

if ~exist([Path '/.fig Files'],'dir'), mkdir([Path '/.fig Files']); end

if ~SRI
    DIR = [Path '/.mat Files/' strjoin(FileName,'') '_' MID '_FULL.mat'];
else
    DIR = [Path '/.mat Files/' strjoin(FileName,'') '_' MID '_SRI.mat'];
end
end