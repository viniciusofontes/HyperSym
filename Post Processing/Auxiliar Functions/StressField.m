% StressField -------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function StressField(MATDIR,FileName,MID,SRI,XYZ,LE,Type,Coord)
% Plot the stresses according to "Type".

load(MATDIR,'DISP','CAUCHY');

persistent NodalCAUCHY

CAUCHY = StressRecovery(CAUCHY(:,:,end),XYZ,LE);

if strcmp(Coord,'CYL') && contains(FileName,'Cylinder')
    if isempty(NodalCAUCHY), NodalCAUCHY = Cart2Cyl(CAUCHY,XYZ); end
    
    % ---------------------------------------------------------------------
    % Correction to deal with boundary effect (at the inner radius near the
    % supports).
    nd = sqrt(size(LE,1));
    
    NodalCAUCHY(1:3,1:(nd+1)) = ...
        NodalCAUCHY(1:3,(nd+2):2*(nd+1));
    NodalCAUCHY(1:3,(nd+1)^2-nd:(nd+1)^2) = ...
        NodalCAUCHY(1:3,(nd+2):2*(nd+1));
    NodalCAUCHY(1:3,(nd+1)^2+1:(nd+1)^2+(nd+1)) = ...
        NodalCAUCHY(1:3,(nd+2):2*(nd+1));
    NodalCAUCHY(1:3,end-nd:end) = ...
        NodalCAUCHY(1:3,(nd+2):2*(nd+1));
    % ---------------------------------------------------------------------
    
    save(MATDIR,'NodalCAUCHY','-append');
elseif strcmp(Coord,'CART')
    NodalCAUCHY = CAUCHY;
else
    error('Wrong coordinate system.');
end

if strcmp(Type,'1')
    row = 1;       StressComp = '$\sigma_1$';
    S = NodalCAUCHY(row,:);
    StressName = 'S1';
elseif strcmp(Type,'2')
    row = 2;       StressComp = '$\sigma_2$';
    S = NodalCAUCHY(row,:);
    StressName = 'S2';
elseif strcmp(Type,'3')
    row = 3;       StressComp = '$\sigma_3$';
    S = NodalCAUCHY(row,:);
    StressName = 'S3';
elseif strcmp(Type,'12') || strcmp(Type,'21')
    row = 4;       StressComp = '$\sigma_{12}$';
    S = NodalCAUCHY(row,:);
    StressName = 'S12';
elseif strcmp(Type,'13') || strcmp(Type,'31')
    row = 5;       StressComp = '$\sigma_{13}$';
    S = NodalCAUCHY(row,:);
    StressName = 'S13';
elseif strcmp(Type,'23') || strcmp(Type,'32')
    row = 6;       StressComp = '$\sigma_{23}$';
    S = NodalCAUCHY(row,:);
    StressName = 'S23';
elseif strcmp(Type,'Mises')
    StressComp = '$\sigma_{vM}$';
    TmpSigma = NodalCAUCHY;
    
    S = sqrt(1/2*((TmpSigma(1,:) - TmpSigma(2,:)).^2 + ...
        (TmpSigma(2,:) - TmpSigma(3,:)).^2 + ...
        (TmpSigma(1,:) - TmpSigma(3,:)).^2 + ...
        6*TmpSigma(4,:).^2 + 6*TmpSigma(5,:).^2 + 6*TmpSigma(6,:).^2)); 
    StressName = 'SvM';
elseif strcmp(Type,'Tresca')
    StressComp = '$\sigma_{Tr}$';
    TmpSigma = NodalCAUCHY;
    
    M = [TmpSigma(1,:) TmpSigma(4,:) TmpSigma(6,:);
         TmpSigma(4,:) TmpSigma(2,:) TmpSigma(5,:);
         TmpSigma(6,:) TmpSigma(5,:) TmpSigma(3,:)];
     
    S = zeros(size(XYZ,1),1);
    
    for n = 1:size(XYZ,1)
        A = M(3*(n-1)+1:3*n,:);
        S(n) = real((max(eig(A)) - min(eig(A)))/2);
    end
    
     StressName = 'STr';
else
    error('Invalid stress field.')
end

PlotStresses(XYZ,LE,DISP,S)
SetColorbar(StressComp,S)

% Save plot(s)
Folder = replace(FileName{1},'_',' ');
subfolder = sprintf('/Output/%s/.fig Files/Stress Distribution/',Folder);

if ~SRI
    IntName = [MID '_FULL'];
else
    IntName = [MID '_SRI'];
end

if length(FileName) == 1
    SaveFIG(gcf,[StressName '_' IntName],subfolder,'.fig')
else
    SaveFIG(gcf,[StressName '_' FileName{2} IntName],subfolder,'.fig')
end
end
%% PLOTSTRESSES
function PlotStresses(Nodes,Elem,DISP,Sigma)

figure; axis off; axis tight; axis equal;  view([0 45]); colormap jet
      
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

UU = DISP(:,:,end);

patch('Faces',cell2mat(Faces),...
    'Vertices',Nodes + UU,...
    'CData',Sigma,...
    'FaceColor','interp',...
    'FaceAlpha',0.5,...
    'EdgeColor','k',...
    'LineStyle','-',...
    'LineWidth',0.7);
pause(1e-06)
end
%% COLORBAR
function SetColorbar(StressComp,Sigma)
% SetColorbar sets the colorbar.

cbar = colorbar;
num_pts = 10;                % Number of points dividing the colorbar

% Colorbar data:
cbar.Label.FontSize = 11;
cbar.Label.FontName = 'Times New Roman';
cbar.Label.HorizontalAlignment = 'center';
cbar.XLabel.VerticalAlignment = 'middle';
cbar.Label.Interpreter = 'latex';
cbar.Limits = [min(Sigma) max(Sigma)];
cbar.Location = 'eastoutside';
cbar.XLabel.String = [StressComp ' (MPa)'];
cbar.XLabel.Rotation = 90;
cbar.XLabel.Position = [3 mean(cbar.Limits) 0];

limits = get(cbar,'Limits');
div = linspace(limits(1),limits(2),num_pts);
cbar.Ticks = div;

if abs(div(end) - div(1)) >= 10
    div = arrayfun(@(x) sprintf('%3.0f',x),div,'un',0);
else
    div = arrayfun(@(x) sprintf('%3.2f',x),div,'un',0);
end

cbar.TickLabels = div;
end