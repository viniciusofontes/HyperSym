% PlotMesh_Strip ----------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x---------------
function PlotMesh_Strip(FileName,Nodes,Elem,Supp,Load)
% Plot mesh in Lagrangian (undeformed) configuration.

figure; grid on; axis off; axis tight; axis equal; hold on; view([55 25])
      
Faces = cell(size(Elem,1),1);

for el = 1:size(Elem,1)
    a = Elem(el,1); b = Elem(el,2); c = Elem(el,3); d = Elem(el,4);
    e = Elem(el,5); f = Elem(el,6); g = Elem(el,7); h = Elem(el,8);
    
    Faces{el} = [a b c d;
                 e f g h;
                 a b f e;
                 d c g h;
                 b c g f;
                 a d h e];
end
      
h = patch('Faces',cell2mat(Faces),...
    'Vertices',Nodes,...
    'FaceColor','white',...
    'FaceAlpha',0.5,...
    'EdgeColor','k',...
    'LineStyle','-',...
    'LineWidth',0.7);
pause(1e-06)

PreDisp  = Supp(Supp(:,3) ~= 0,:,:);   % non-zero displacements
NullDisp = Supp(Supp(:,3) == 0,:,:);   % zero displacements

h1 = plot3(Nodes(NullDisp(:,1),1),Nodes(NullDisp(:,1),2),Nodes(NullDisp(:,1),3),'r>',...
    'MarkerSize',5,...
    'MarkerFaceColor','red');
Dofs = 'Fixed Dofs';

if ~isempty(PreDisp)
    h2 = plot3(Nodes(PreDisp(:,1),1),Nodes(PreDisp(:,1),2),Nodes(PreDisp(:,1),3),'g>',...
        'MarkerSize',5,...
        'MarkerFaceColor','green');
    NonzeroBdCnd = 'Free Dofs';
end

if ~isempty(Load)
    h2 = plot3(Nodes(Load(:,1),1),Nodes(Load(:,1),2),Nodes(Load(:,1),3),'g>',...
        'MarkerSize',5,...
        'MarkerFaceColor','green');
    NonzeroBdCnd = 'Apllied Load';
end

legend([h1,h2],Dofs,NonzeroBdCnd,'Location','nw')

SaveMesh(h,FileName,'.fig')
end
%% SAVE MESH
function SaveMesh(h,FileName,ext)
% SaveMesh saves the mesh in the current working directory.

Folder = replace(FileName{1},'_',' '); % Eliminate the part after '_' 
Path = [pwd '/Output/' Folder '/.fig Files/'];   % Define output folder

if exist(Path,'dir') == 0, mkdir(Path); end

if length(FileName) == 1
    saveas(h,[Path 'Mesh&BC' ext]);
else
    saveas(h,[Path ['Mesh&BC' FileName{2}] ext]);
end

fprintf('Mesh figure saved!\n\n')
end