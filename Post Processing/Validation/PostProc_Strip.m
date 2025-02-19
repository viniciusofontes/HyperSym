% PostProc_Strip ----------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function PostProc_Strip(DIR,SRI,XYZ,LE,BdBox,N)
% Postprocess data in the strip example.

if SRI, IntName = 'SRI'; else, IntName = 'FULL'; end

NSTEPAnsys = 20;   % number of time steps used in ANSYS

MAT = 'Mooney-Rivlin';
%% Ansys data
% Convergence Data --------------------------------------------------------
lineDispI = 8;     % initial line where the displacement is stored
lineStress = 27;   % initial line where the stress is stored
DispI = zeros(length(N),3);
CauchyI = zeros(length(N),8);

Ansys.IntForce = zeros(length(N),3);

Ansysfolder = [pwd '/Ansys/Strip/Output'];

for a = 1:length(N)
    % Displacement and stress in reference nodes
    FileNameI = sprintf(...
        'Results_Selected_Nodes_%s_Mesh=%dx%dx1.txt',MAT,N(a),N(a));   
    fullFileNameI = fullfile(Ansysfolder,FileNameI);
    
    [CauchyI(a,:),DispI(a,:)] = ReadAnsysDataI(fullFileNameI,lineStress,lineDispI);
    
    % Internal force
    FileNameII = sprintf(...
        'Results_IntForce_%s_Mesh=%dx%dx1.txt',MAT,N(a),N(a));
    fullFileNameII = fullfile(Ansysfolder,FileNameII);
    
    Ansys.IntForce(a,:) = ReadAnsysDataII(fullFileNameII,2);
end

Ansys.U = DispI(DispI(:,2) ~= 0,2);
Ansys.Pos = BdBox(2) + Ansys.U;

Ansys.CauchyMises = CauchyI(:,end-1);

% Equilibrium Path --------------------------------------------------------
lineDispIII = [5 31];
DispIII = zeros(NSTEPAnsys+1,2*3);

for b = 1:2
    FileNameIII = sprintf('Results_Eq_Path_%s.txt',MAT);
    fullFileNameIII = fullfile(Ansysfolder,FileNameIII);
    
    [lbd,DispIII(:,3*b-2:3*b)] = ReadAnsysDataIII(fullFileNameIII,lineDispIII(b));
end

Ansys.LF = lbd;
Ansys.UX = DispIII(:,4);
Ansys.UY = DispIII(:,2);
%% Matlab data
% Convergence Data --------------------------------------------------------
qtd = 2;      zElem = 1;
Matlab.U = zeros(qtd*length(N),3);
Matlab.CauchyMises = zeros(length(N),qtd);
Matlab.IntForce = zeros(length(N),3);

for d = 1:length(N)
    % Nodal stresses in the last time step
    load(DIR{d},'DISP','CAUCHY','Fi','FAC')
    
    U = DISP(:,:,end);       % Nodal displacement in the last time step
        
    Cauchy = StressRecovery(CAUCHY(:,:,end),XYZ{d},LE{d});
    
    CauchyMises = sqrt(1/2*((Cauchy(1,:) - Cauchy(2,:)).^2 + ...
        (Cauchy(2,:) - Cauchy(3,:)).^2 + ...
        (Cauchy(1,:) - Cauchy(3,:)).^2 + ...
        6*Cauchy(4,:).^2 + 6*Cauchy(5,:).^2 + 6*Cauchy(6,:).^2));
    
    nd = [(N(d)+1)^2*zElem - N(d) (N(d)+1)^2*zElem];
    
    Matlab.U(qtd*(d-1)+1:qtd*d,:) = U(nd,:);
    Matlab.CauchyMises(d,:) = CauchyMises(:,nd);
    
    % Internal force on the edge elements
    RightNodes = abs(XYZ{d}(:,1) - BdBox(2)) < eps;
    
    EdgeIntForce = sum(Fi(RightNodes,:,end),1);
    Matlab.IntForce(d,:) = EdgeIntForce;
end

Matlab.U = Matlab.U(1:qtd:end,2);
Matlab.Pos = BdBox(2) + Matlab.U;

% Equilibrium Path --------------------------------------------------------
BaseLeftTopNode = abs(XYZ{end}(:,1) - 0) < eps &...
    abs(XYZ{end}(:,2) - BdBox(2)) < eps & ...
    abs(XYZ{end}(:,3) - 0) < eps;

BaseRightBottomNode = abs(XYZ{end}(:,1) - BdBox(2)) < eps &...
    abs(XYZ{end}(:,2) - 0) < eps & ...
    abs(XYZ{end}(:,3) - 0) < eps;

UX = DISP(BaseRightBottomNode,1,:);
UY = DISP(BaseLeftTopNode,2,:);
Matlab.UY = [0 reshape(UY,[],length(FAC))];
Matlab.UX = [0 reshape(UX,[],length(FAC))];
Matlab.LF = [0 FAC];
%% Comparison
% Gap in Cauchy Mises stress for different meshes -------------------------
MisesGap = Matlab.CauchyMises(:,2) - Ansys.CauchyMises;
fprintf('Gap in Mises data:\n'); disp(mean(MisesGap));
%% Plot
% Convergence data --------------------------------------------------------
h1 = figure; hold on; box on;

yyaxis left
axl = gca;
axl.YColor = 'b';

plot(log2(N.^2),Matlab.Pos,...
    'LineStyle','-',...
    'LineWidth',1,...
    'Marker','s',...
    'MarkerSize',6,...
    'Color','b',...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0 255 255]/255);

plot(log2(N.^2),Ansys.Pos,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Marker','o',...
    'MarkerSize',6,...
    'Color','b',...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0 255 255]/255);

ylabel('Final Position (mm)')

yT = get(gca,'YTick'); ylim([yT(1) yT(end)])
ytickformat('%.2f');   yticks(linspace(yT(1),yT(end),7))

yyaxis right;
axr = gca;
axr.YColor = [0 176 80]/255;

plot(log2(N.^2),Matlab.IntForce(:,1),...
    'LineStyle','-',...
    'LineWidth',1,...
    'Marker','s',...
    'MarkerSize',6,...
    'Color',[0 176 80]/255,...
    'MarkerEdgeColor',[0 176 80]/255,...
    'MarkerFaceColor','g')

plot(log2(N.^2),Ansys.IntForce(:,1),...
    'LineStyle','--',...
    'LineWidth',1,...
    'Marker','o',...
    'MarkerSize',6,...
    'Color',[0 176 80]/255,...
    'MarkerEdgeColor',[0 176 80]/255,...
    'MarkerFaceColor','g')

xlabel('Number of Elements')
ylabel('Edge Force on the Right Face (N)')

yT = get(gca,'YTick'); ylim([yT(1) yT(end)])
ytickformat('%.2f');   yticks(linspace(yT(1),yT(end),7))

xticks(log2(N.^2));     xT = get(gca,'XTick');
xticklabels(cellstr(num2str(xT(:),'$2^{%d}$')))
xlim([log2(N(1)^2)-0.1 log2(N(end)^2)+0.1])

if SRI, legend('HyperSym','ANSYS','HyperSym','ANSYS','Location','northeast')
else,   legend('HyperSym','ANSYS','HyperSym','ANSYS','Location','east')
end

% Save plot(s)
GraphName = 'Convergence_Data'; 
subfolder = '/Output/Strip/.fig Files/';
SaveFIG(h1,[GraphName sprintf('_%s_%s',MAT,IntName)],subfolder,'.fig');

% Stresses ----------------------------------------------------------------
h2 = figure; hold on; box on

plot(log2(N.^2),Matlab.CauchyMises(:,2),...
    'LineStyle','-',...
    'LineWidth',1,...
    'Marker','s',...
    'MarkerSize',6,...
    'Color',[220 20 60]/255,...
    'MarkerEdgeColor',[220 20 60]/255,...
    'MarkerFaceColor',[255 192 203]/255)

plot(log2(N.^2),Ansys.CauchyMises,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Marker','o',...
    'MarkerSize',6,...
    'Color',[220 20 60]/255,...
    'MarkerEdgeColor',[220 20 60]/255,...
    'MarkerFaceColor',[255 192 203]/255)

xlabel('Number of Elements')
ylabel('$\sigma_{vM}$ (MPa)')

xticks(log2(N.^2));     xT = get(gca,'XTick');
xticklabels(cellstr(num2str(xT(:),'$2^{%d}$')))
xlim([log2(N(1)^2)-0.1 log2(N(end)^2)+0.1])

legend('HyperSym','ANSYS','Location','southeast')

% Save plot(s)
GraphName = 'Stresses'; 
subfolder = '/Output/Strip/.fig Files/';
SaveFIG(h2,[GraphName sprintf('_%s_%s',MAT,IntName)],subfolder,'.fig');

% Equilibrium Path --------------------------------------------------------
h3 = figure; hold on; box on
axis([-35 205 0 1])

plot(Matlab.UX,Matlab.LF,...
    'LineStyle','-',...
    'Color',[220 20 60]/255,...
    'LineWidth',1)
plot(Matlab.UY,Matlab.LF,...
    'LineStyle','-',...
    'Color',[0 0 1],...
    'LineWidth',1)
plot(Ansys.UX,Ansys.LF,...
    'LineStyle','none',...
    'LineWidth',1,...
    'Marker','o',...
    'MarkerSize',6,...
    'Color',[220 20 60]/255,...
    'MarkerEdgeColor',[220 20 60]/255,...
    'MarkerFaceColor',[255 192 203]/255)
plot(Ansys.UY,Ansys.LF,...
    'LineStyle','none',...
    'LineWidth',1,...
    'Marker','o',...
    'MarkerSize',6,...
    'Color','b',...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0 255 255]/255)

xlabel('$u_i$ (mm)')    % displacement
ylabel('Load Factor')   % load factor
legend('$u_1$, HyperSym','$u_2$, HyperSym',...
    '$u_1$, ANSYS','$u_2$, ANSYS',...
    'Location','southeast')

annotation(gcf,'textarrow',[0.233 0.180],[0.543 0.495],...
    'String',{'A'},...
    'Color',[0 0 1],...
    'FontSize',10,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','baseline',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','plain',...
    'HeadWidth',6,...
    'HeadLength',7);

annotation(gcf,'textarrow',[0.499 0.552],[0.549 0.505],...
    'String',{'B'},...
    'Color',[220 20 60]/255,...
    'FontSize',10,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','baseline',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','plain',...
    'HeadWidth',6,...
    'HeadLength',7);

xticks([-32,-16,0,17:17:204]);

% Save plot(s)
GraphName = 'Equilibrium_Path'; 
subfolder = '/Output/Strip/.fig Files/';
SaveFIG(h3,[GraphName sprintf('_%s_%s',MAT,IntName)],subfolder,'.fig');
end
%% ReadAnsysData I
function [CauchyStress,varargout] = ReadAnsysDataI(fileName,lineStress,lineDisp)
% ReadAnsysDataI reads an output file generated by Ansys.

StressID = fopen(fileName,'r');

item = textscan(StressID,'%f %f %f %f %f %f %f %f',...
    'Delimiter','\n','Headerlines',lineStress-1);
CauchyStress = cell2mat(item);
fclose(StressID);

n = nargout;

if n == 2
    DispID = fopen(fileName,'r');
    
    item = textscan(DispID,'%f %f %f',1,...
        'Delimiter','\n','Headerlines',lineDisp-1);
    varargout{1} =  cell2mat(item);
    fclose(DispID);
end
end
%% ReadAnsysData II
function Fi = ReadAnsysDataII(fileName,lineFi)
% ReadAnsysDataII reads an output file generated by Ansys.

DispID = fopen(fileName,'r');

item = textscan(DispID,'%f %f %f %*f %*f %*f',...
    'Delimiter','\n','Headerlines',lineFi-1);
Fi = cell2mat(item);
fclose(DispID);
end
%% ReadAnsysData III
function [lbd,Disp] = ReadAnsysDataIII(fileName,lineDisp)
% ReadAnsysDataIII reads an output file generated by Ansys.

DispID = fopen(fileName,'r');

item = textscan(DispID,'%f %f %f %f',...
    'Delimiter','\n','Headerlines',lineDisp-1);
lbd = item{1};
Disp = [item{:,2} item{:,3} item{:,4}];
fclose(DispID);
end