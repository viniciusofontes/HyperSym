% PostProc_UniaxialTraction_Comp_Models -----------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function PostProc_UniaxialTraction_Comp_Models(DIR,SRI,XYZ,LE,NREF,MID,E,nu)
% Postprocess data in the uniaxial traction made of compressible models
% example.

[h1,h2] = PlotAnalyticalStresses(E,nu,MID);

if SRI, IntName = 'SRI'; else, IntName = 'FULL'; end

NSTEP = 20;
%% Matlab data
Marker = {'d','o','^','x','o','^','x'};
Color = {[1 0.8 0]; [0 0 1]; [0 0 1]; [0 0 1];
    [1 0 0]; [1 0 0]; [1 0 0];[0.5 0.5 0.5]};
LineStyle = {'none','none','none','none','none','none','none'};
MarkerSize = 6*ones(7,1);
%% Plot
% Stress vs. stretch ------------------------------------------------------
U = zeros(NSTEP+1,3,length(MID));
Matlab.Cauchy = zeros(NSTEP+1,6,length(MID));

figure(h1); hold all; box on; axis tight

for m = 1:4
    load(DIR{m},'DISP','CAUCHY');

    for t = 1:NSTEP
        U(t+1,:,m) = DISP(NREF,:,t);  % nodal dis in the last time step
        
        [~,S] = StressRecovery(CAUCHY(:,:,t),XYZ,LE);
        Matlab.Cauchy(t+1,:,m) = S;
    end  
    
    Matlab.StretchX = 1 + U(:,1,m);
    
    plot(Matlab.StretchX,Matlab.Cauchy(:,1,m)/1E0,...
        'Color',Color{m},...
        'LineStyle',LineStyle{m},...
        'LineWidth',1,...
        'Marker',Marker{m},...
        'MarkerSize',MarkerSize(m),...
        'DisplayName',['HyperSym, ' MID{m}])
end

xlabel('$\gamma$')            % Stretch ratio on x-direction
ylabel('$\sigma_{11}$ (MPa)') % Cauchy stress

% Joining analytical and Matlab datas -------------------------------------
MinRange = 0.9;                    MaxRange = 1.1;

axObj1  = get(h1,'Children');
datObj1 = get(axObj1,'Children');

x1 = get(datObj1{2},'XData');
y1 = get(datObj1{2},'YData');
Color1  = get(datObj1{2},'Color');
LineStyle1  = get(datObj1{2},'LineStyle');
Marker1 = get(datObj1{2},'Marker');
MarkerTam1  = get(datObj1{2},'MarkerSize');

% Zoom Image
annotation(h1,'rectangle',...
    [0.432104166666667 0.53023758099352 0.0850833333333334 0.0453563714902808]);
annotation(h1,'line',...
    [0.431770833333333 0.1802],...
    [0.531317494600432 0.6112]);
annotation(h1,'line',...
    [0.516666666666667 0.4052],...
    [0.575593952483801 0.8747]);

Pos = get(gca,'Position');
ax = axes('Position',[Pos(1)+0.05 Pos(2)+0.5 Pos(3)-0.55 Pos(4)-0.55]);

FigInFig(ax,[MinRange MaxRange;-200 200],x1,y1,...
    Color1,LineStyle1,Marker1,MarkerTam1)
% -------------------------------------------------------------------------

% Save plot(s)
GraphName = 'Stress_x_Stretch_SVK';
subfolder = '/Output/Uniaxial Traction Comp Models/.fig Files/';
SaveFIG(gcf,[GraphName sprintf('_%s',IntName)],subfolder,'.fig')

figure(h2); hold all; box on; axis tight

for m = 5:length(MID)
    load(DIR{m},'DISP','CAUCHY');

    for t = 1:NSTEP
        U(t+1,:,m) = DISP(NREF,:,t);  % nodal dis in the last time step
        
        [~,S] = StressRecovery(CAUCHY(:,:,t),XYZ,LE);
        Matlab.Cauchy(t+1,:,m) = S;
    end  
    
    Matlab.StretchX = 1 + U(:,1,m);
    
    plot(Matlab.StretchX,Matlab.Cauchy(:,1,m)/1E0,...
        'Color',Color{m},...
        'LineStyle',LineStyle{m},...
        'LineWidth',1,...
        'Marker',Marker{m},...
        'MarkerSize',MarkerSize(m),...
        'DisplayName',['HyperSym, ' MID{m}])
end

xlabel('$\gamma$')            % Stretch ratio on x-direction
ylabel('$\sigma_{11}$ (MPa)') % Cauchy stress

% Joining analytical and Matlab datas -------------------------------------
% MinRange = 0.9;                    MaxRange = 1.1;

axObj2  = get(h2,'Children');
datObj2 = get(axObj2,'Children');

x2 = get(datObj2{2},'XData'); 
y2 = get(datObj2{2},'YData');
Color2  = get(datObj2{2},'Color');
LineStyle2  = get(datObj2{2},'LineStyle');
Marker2 = get(datObj2{2},'Marker');
MarkerTam2  = get(datObj2{2},'MarkerSize');

% Zoom Image
annotation(h2,'rectangle',...
    [0.431770833333333 0.793736501079913 0.0859375 0.0362634989200858]);
annotation(h2,'line',...
    [0.431770833333333 0.605729166666667],...
    [0.795896328293736 0.48488120950324]);
annotation(h2,'line',...
    [0.517708333333333 0.8296875],...
    [0.830453563714903 0.750539956803456]);

Pos = get(gca,'Position');
ax = axes('Position',[Pos(1)+0.475 Pos(2)+0.375 Pos(3)-0.55 Pos(4)-0.55]);

FigInFig(ax,[MinRange MaxRange;-200 200],x2,y2,...
    Color2,LineStyle2,Marker2,MarkerTam2)
% -------------------------------------------------------------------------

% Save plot(s)
GraphName = 'Stress_x_Stretch_nH';
subfolder = '/Output/Uniaxial Traction Comp Models/.fig Files/';
SaveFIG(gcf,[GraphName sprintf('_%s',IntName)],subfolder,'.fig')
end