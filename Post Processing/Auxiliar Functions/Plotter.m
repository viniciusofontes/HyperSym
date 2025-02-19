% Plotter -----------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function Plotter(XYZ,LE,DIR,ARGS)
% Plot mesh.

Faces = cell(size(LE,1),1);

for IE = 1:size(LE,1)
    a = LE(IE,1); b = LE(IE,2); c = LE(IE,3); d = LE(IE,4);
    e = LE(IE,5); f = LE(IE,6); g = LE(IE,7); h = LE(IE,8);
    
    Faces{IE} = [...
        a b c d;
        e f g h;
        a b f e;
        d c g h;
        a d h e;
        b c g f];
end

figure;
load(DIR,'DISP');

for ITER = 1:size(DISP,3)   
    clf; hold on; grid on; box on; axis equal; 
    axis(ARGS.axis); view(ARGS.view)
    xlabel('$x_1$'); ylabel('$x_2$'); zlabel('$x_3$','Rotation',0);
    
    patch('Faces',cell2mat(Faces),...
        'Vertices',XYZ,...
        'FaceColor','none',...
        'EdgeColor','k',...
        'LineStyle','--',...
        'LineWidth',1);
    
    patch('Faces',cell2mat(Faces),...
        'Vertices',XYZ + DISP(:,:,ITER),...
        'FaceColor','blue',...
        'FaceAlpha',0.5,...
        'EdgeColor','k',...
        'LineStyle','-',...
        'LineWidth',1);
    
    if isfield(ARGS,'xticks'), xticks(ARGS.xticks); end
    if isfield(ARGS,'yticks'), yticks(ARGS.yticks); end
    if isfield(ARGS,'zticks'), zticks(ARGS.zticks); end
    pause(5E-2)
end
end