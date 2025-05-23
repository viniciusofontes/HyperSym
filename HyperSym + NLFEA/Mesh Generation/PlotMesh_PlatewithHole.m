% PlotMesh_PlatewithHole --------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function PlotMesh_PlatewithHole(FileName,Nodes,Elem,Supp,Load)
% Plot mesh in Lagrangian (undeformed) configuration.

figure; grid on; axis off; axis tight; axis equal; hold on; view([-20 40])
      
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

legend([h1,h2],Dofs,NonzeroBdCnd,'Location','se')

SaveMesh(h,FileName,'.fig')
end
%% SAVE MESH
function SaveMesh(h,FileName,ext)
% SaveMesh saves the mesh in the current working directory.

Folder = replace(FileName,'_',' ');    % Save output on a folder without '_' 
Path = [pwd '/Output/' Folder '/.fig Files/'];   % Defining output folder

if exist(Path,'dir') == 0, mkdir(Path); end
                 
saveas(h,[Path 'Mesh&BC' ext]);

fprintf('Mesh figure saved!\n\n')
end