% SaveFig -----------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function SaveFIG(h,FileName,subfolder,ext)
% Adjust to save a figure as a '.pdf' file.
% The first commands adjusts the figure size.
% The following ones save the figure in the desired format given by 'ext'.
% The last part checks whether the folder exists or not and move the file
% to the correct folder.

if ~strcmp(ext,'.fig') || ~strcmp(ext,'fig'), ext = '.fig'; end

set(h,'Units','inches');
screenposition = get(h,'Position');
set(h,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

if exist([pwd subfolder],'dir') == 0, mkdir([pwd '\' subfolder]); end

saveas(h,[pwd subfolder FileName ext]);

fprintf('Figure %i saved!\n',h.Number)
end