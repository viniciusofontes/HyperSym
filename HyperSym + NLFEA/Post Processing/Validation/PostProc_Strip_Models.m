% PostProc_Strip_Models ---------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function PostProc_Strip_Models(DIR,SRI,XYZ,LE,BdBox,MID,State)
% Postprocess data in the strip_models example.

if SRI, IntName = 'SRI'; else, IntName = 'FULL'; end

XYZ_max = max(XYZ);
XYZ_min = min(XYZ);
nodeA = find(XYZ(:,1) == XYZ_min(1));
nodeA = nodeA(XYZ(nodeA,2) == XYZ_max(2));
nodeA = nodeA(XYZ(nodeA,3) == XYZ_min(3));
%% Matlab data
% Equilibrium Path --------------------------------------------------------
[Matlab.U,Matlab.LF,Matlab.VOL] = deal(cell(length(MID),1));

for m = 1:length(MID)
    load(DIR{m},'DISP','FAC')
        
    Matlab.U{m} = [0; squeeze(DISP(nodeA,2,:))]; % Nodal displacement
    Matlab.LF{m} = [0; FAC'];                    % Load factor
    
    Matlab.VOL{m} = ComputeVolume(XYZ,LE,DISP);
end
%% Plot
% Equilibrium Path --------------------------------------------------------
h1 = figure; hold on; box on; axis tight

p.Marker = {'d','o','^','x','o','^','x','s'};
p.Color = {[1 0.8 0]; [0 0 1]; [0 0 1]; [0 0 1];
    [1 0 0]; [1 0 0]; [1 0 0]; [0 0.8 0]};

for m = 1:length(MID)
    plot(Matlab.U{m},Matlab.LF{m},...
        'LineStyle','-',...
        'LineWidth',1,...
        'Marker',p.Marker{m},...
        'MarkerSize',6,...
        'Color',p.Color{m},...
        'DisplayName',MID{m});
end

xlabel('$u_2$ (mm)');   ylabel('Load Factor')

ax = gca;
xlim([ax.XLim(1) 0]);   ylim([0 ax.YLim(end)])

legend('-DynamicLegend','Location','southwest')

% Save plot(s)
GraphName = 'Equilibrium_Path';
subfolder = '/Output/Strip Models/.fig Files/';
SaveFIG(h1,[GraphName sprintf('_%s',IntName)],subfolder,'.fig')

% Volume Ratio ------------------------------------------------------------
h2 = figure; hold on; box on; axis tight

plot([0 1],[1 1]*100,...
    'LineStyle','--',...
    'LineWidth',1,...
    'Color','k');

annotation(h2,'textarrow',[0.480 0.444],[0.567 0.689],...
    'String',{'Original Volume'},...
    'Color',[0 0 0],...
    'FontSize',10,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','plain',...
    'HeadWidth',6,...
    'HeadLength',7);

H = zeros(length(MID),1);

for m = 1:length(MID)
    H(m) = plot(Matlab.LF{m},Matlab.VOL{m}/Matlab.VOL{m}(1)*100,...
        'LineStyle','-',...
        'LineWidth',1,...
        'Marker',p.Marker{m},...
        'MarkerSize',6,...
        'Color',p.Color{m},...
        'DisplayName',MID{m});
end

xlabel('Load Factor');  ylabel('Volume Ratio (\%)')

ax = gca;
xlim([0 ax.XLim(end)]); ylim([0 round(ax.YLim(end),2,'significant')])

legend(H,'Position',[0.646 0.236 0.241 0.328])

% Save plot(s)
GraphName = 'Volume';
subfolder = '/Output/Strip Models/.fig Files/';
SaveFIG(h2,[GraphName sprintf('_%s',IntName)],subfolder,'.fig')

% Deformed Configuration --------------------------------------------------
DeformedConfig(DIR,IntName,XYZ,BdBox,MID,p,State)
end
%% Volume
function VOL = ComputeVolume(XYZ,LE,DISP)
% ComputeVolume computes the volume at every time step for all elements in
% the FE model.

N = size(DISP,3) + 1;   VOL = zeros(N,1);   NE = size(LE,1);

for n = 1:N
    VOLUME = 0;

    for I = 1:NE
        if n == 1
            ELXY = XYZ(LE(I,:),:);
        else
            ELXY = XYZ(LE(I,:),:) + DISP(LE(I,:),:,n-1);
        end
        
        DET = Det([0 0 0],ELXY);
        DVOL = 8*DET;
        
        VOLUME = VOLUME + DVOL;
    end
    
    VOL(n) = VOLUME;
end
end
%% Determinat
function DET = Det(XI,ELXY)
% Compute determinant of hexahedron element.

XNODE = [...
    -1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
QUAR = 0.125;

SF = zeros(8,1);        DSF = zeros(3,8);

for I = 1:8
    XP = XNODE(1,I);
    YP = XNODE(2,I);
    ZP = XNODE(3,I);
    
    XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];
    
    SF(I) = QUAR*XI0(1)*XI0(2)*XI0(3);
    DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
    DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
    DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
end

GJ  = DSF*ELXY; 
DET = det(GJ);
end