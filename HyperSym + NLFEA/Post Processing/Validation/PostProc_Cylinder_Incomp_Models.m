% PostProc_Cylinder_Incomp_Models -----------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function PostProc_Cylinder_Incomp_Models(DIR,SRI,XYZ,LE,BdBox,MID)
% Postprocess data in the thich-walled cylinder made of incompressible 
% models example.

if SRI, IntName = 'SRI'; else, IntName = 'FULL'; end

MAT = {'Mooney-Rivlin'; 'Yeoh'};
%% Ansys data
Ansysfolder = [pwd '/Ansys/Thick Walled Cylinder/Output'];
CS = 'CYL';

% Equilibrium Path --------------------------------------------------------
lineDisp = [5; 21];
Disp = zeros(11,6);
Ansys.CylCoord.Path.AvgCauchy = cell(length(MID),1);
Ansys.CylCoord.Path.Position = cell(length(MID),1);

for m = 1:length(MID)
    FileName1 = sprintf('Results_Eq_Path_%s_%s.txt',MAT{m},CS);
    fullFileName1 = fullfile(Ansysfolder,FileName1);
    
    for i = 1:2
        [lbd,Disp(:,3*i-2:3*i)] = ReadAnsysDataII(fullFileName1,lineDisp(i));
    end
    
    % Average of nodal displacement for the bottom and top data
    AvgU = [...
        mean([Disp(:,1) Disp(:,4)],2) ...
        mean([Disp(:,2) Disp(:,5)],2) ...
        mean([Disp(:,3) Disp(:,6)],2)];
    
    Ansys.LF(:,m) = lbd;
    Ansys.CylCoord.AvgU{m} = AvgU;
    
    % Stress Distribution -------------------------------------------------
    lineStress = 3;
    
    FileName2 = sprintf('Results_BottomStress_%s_%s.txt',MAT{m},CS);
    fullFileName2 = fullfile(Ansysfolder,FileName2);
    
    FileName3 = sprintf('Results_TopStress_%s_%s.txt',MAT{m},CS);
    fullFileName3 = fullfile(Ansysfolder,FileName3);
    
    % Bottom Path
    Ansys.CylCoord.BottomPath.Cauchy = ReadAnsysDataI(fullFileName2,lineStress);
    Ansys.CylCoord.BottomPath.Position = ...
        (linspace(BdBox(1),BdBox(2),size(Ansys.CylCoord.BottomPath.Cauchy,1)))';
    
    % Top Path
    Ansys.CylCoord.TopPath.Cauchy = ReadAnsysDataI(fullFileName3,lineStress);
    Ansys.CylCoord.TopPath.Position = Ansys.CylCoord.BottomPath.Position;
    
    % Average of nodal stress for the bottom and top data
    for i = 1:size(Ansys.CylCoord.BottomPath.Cauchy,1)
        Ansys.CylCoord.Path.AvgCauchy{m}(i,:) = ...
            mean([Ansys.CylCoord.BottomPath.Cauchy(i,:); ...
            Ansys.CylCoord.TopPath.Cauchy(i,:)]);
    end

    Ansys.CylCoord.Path.Position{m} = Ansys.CylCoord.BottomPath.Position';
end
%% Matlab data
% Equilibrium Path --------------------------------------------------------
BottomNode = find(XYZ(:,1) == XYZ(:,2) &...
     abs(XYZ(:,3) - 0)  < eps);

TopNode = find(XYZ(:,1) == XYZ(:,2) &...
     abs(XYZ(:,3) - BdBox(6))  < eps);
 
Matlab.CylCoord.AvgU = cell(length(MID),1);
Matlab.CylCoord.Path.AvgCauchy = cell(length(MID),1);
Matlab.CylCoord.Path.Position = cell(length(MID),1);
 
 for m = 1:length(MID)
     load(DIR{m},'DISP','CAUCHY','FAC')
     
     NSTEP = size(DISP,3);

     UBottom = DISP(BottomNode,:,:);
     URefBottom = UBottom(1,:,:);
     URefBottom = reshape(URefBottom,[],length(FAC));
     UTop = DISP(TopNode,:,:);
     URefTop = UTop(1,:,:);
     URefTop = reshape(URefTop,[],length(FAC));
     
     Matlab.CartCoord.AvgU = zeros(NSTEP+1,3);
     Matlab.CartCoord.Cauchy = zeros(6,length(XYZ),NSTEP+1);
     
     for t = 1:NSTEP
         Matlab.CartCoord.AvgU(t+1,:) = mean([URefBottom(:,t) URefTop(:,t)],2);
         
         Matlab.CartCoord.Cauchy(:,:,t+1) = StressRecovery(CAUCHY(:,:,t),XYZ,LE);
     end

     % Matlab.CylCoord.Ur = vecnorm(Matlab.CartCoord.U,2,2);
     [theta,r,z] = cart2pol(Matlab.CartCoord.AvgU(:,1),...
         Matlab.CartCoord.AvgU(:,2),...
         Matlab.CartCoord.AvgU(:,3));
     
     theta(1) = theta(2);  % the first value is set to zero which is wrong,
     ...                     since all values correspond to the same path.
     ...                     The above is a simple corection to it
         
     Matlab.CylCoord.AvgU{m} = [r rad2deg(theta) z];
     Matlab.LF(:,m) = [0 FAC]';
     
     % Stress Distribution ------------------------------------------------
     Matlab.CylCoord.BottomPath.Position = (linspace(BdBox(1),BdBox(2),length(BottomNode)))';
     Matlab.CylCoord.TopPath.Position = Matlab.CylCoord.BottomPath.Position;
     Matlab.CylCoord.Path.Position{m} = Matlab.CylCoord.TopPath.Position;
     
     Matlab.CartCoord.BottomPath.Cauchy = Matlab.CartCoord.Cauchy(:,BottomNode,end);
     Matlab.CartCoord.TopPath.Cauchy = Matlab.CartCoord.Cauchy(:,TopNode,end);
     
     Matlab.CylCoord.BottomPath.Cauchy = ...
         Cart2Cyl(Matlab.CartCoord.BottomPath.Cauchy,XYZ(BottomNode,:));
     Matlab.CylCoord.TopPath.Cauchy = ...
         Cart2Cyl(Matlab.CartCoord.BottomPath.Cauchy,XYZ(TopNode,:));
     
     Matlab.CylCoord.Path.AvgCauchy{m} = mean([Matlab.CylCoord.BottomPath.Cauchy;
         Matlab.CylCoord.TopPath.Cauchy]);
     
%      Matlab.CylCoord.Path.AvgCauchy = ...
%          zeros(size(Matlab.CylCoord.BottomPath.Cauchy));

     for j = 1:size(Matlab.CylCoord.BottomPath.Cauchy,1)
         Matlab.CylCoord.Path.AvgCauchy{m}(j,:) = ...
             mean([Matlab.CylCoord.BottomPath.Cauchy(j,:); ...
             Matlab.CylCoord.TopPath.Cauchy(j,:)]);
     end
 end

% Radial stress at r = 10% of thickness -----------------------------------
fprintf('\n')

for m = 1:length(MID)
    Sr = Matlab.CylCoord.Path.AvgCauchy{m}(1,:)';
    ri = BdBox(1);  ro = BdBox(2);   t = ro - ri;
    r = linspace(ri,ro,length(Sr))'; r_10 = ri + 0.1*t;
    Sr_10 = interp1(r,Sr,r_10);      psi2MPa = 6.8945757293E-3;
    
    fprintf('Model #%g: %s\n',m,MID{m})
    fprintf('----------------------------------------\n\n')
    fprintf('Radial stress at 10%% of the thickness:\n%g MPa\n%g psi\n',...
        Sr_10,Sr_10/psi2MPa)
    fprintf('----------------------------------------\n\n')
end
%% Plot
LineColor = {[0 0.8 0]; [0.5 0.5 0.5]};
EdgeColor = LineColor;
FaceColor = {[0 1 0]; [0.7 0.7 0.7]};

% Equilibrium Path --------------------------------------------------------
h1 = figure; hold on; box on; axis tight
  
for m = 1:length(MID)
    plot(Matlab.CylCoord.AvgU{m}(:,1),Matlab.LF(:,m),...
        'Color',LineColor{m},...
        'LineWidth',1,...
        'DisplayName',['HyperSym, ' MID{m}])
    plot(Ansys.CylCoord.AvgU{m}(:,1),Ansys.LF(:,m),...
        'Color',LineColor{m},...
        'LineStyle','none',...
        'LineWidth',1,...
        'Marker','o',...
        'MarkerSize',6,...
        'MarkerEdgeColor',EdgeColor{m},...
        'DisplayName',['ANSYS, ' MID{m}])
    
    legend('Location','southeast')
end

xlabel('$u_r$ (mm)')    % Displacement
ylabel('Load Factor')   % Load factor

% Save plot(s)
GraphName = 'Equilibrium_Path'; 
subfolder = '/Output/Thick Walled Cylinder Incomp Models/.fig Files/';
SaveFIG(h1,[GraphName sprintf('_%s',IntName)],subfolder,'.fig');

% Stress Distribution -----------------------------------------------------
h2 = cell(length(MID),1);

for m = 1:length(MID)
    h2{m} = figure; hold on; box on; axis tight
    
    plot(Matlab.CylCoord.Path.Position{m},Matlab.CylCoord.Path.AvgCauchy{m}(1,:),...
        'Color',[0 0 1],...
        'LineWidth',1)
    plot(Ansys.CylCoord.Path.Position{m},Ansys.CylCoord.Path.AvgCauchy{m}(:,1),...
        'LineStyle','none',...
        'LineWidth',1,...
        'Marker','o',...
        'MarkerSize',6,...
        'MarkerEdgeColor','b')
    
    plot(Matlab.CylCoord.Path.Position{m},Matlab.CylCoord.Path.AvgCauchy{m}(2,:),...
        'Color',[220 20 60]/255,...
        'LineWidth',1)
    plot(Ansys.CylCoord.Path.Position{m},Ansys.CylCoord.Path.AvgCauchy{m}(:,2),...
        'LineStyle','none',...
        'LineWidth',1,...
        'Marker','^',...
        'MarkerSize',6,...
        'Color',[220 20 60]/255,...
        'MarkerEdgeColor',[220 20 60]/255)
    
    plot(Matlab.CylCoord.Path.Position{m},Matlab.CylCoord.Path.AvgCauchy{m}(3,:),...
        'Color',[60 179 113]/255,...
        'LineWidth',1)
    plot(Ansys.CylCoord.Path.Position{m},Ansys.CylCoord.Path.AvgCauchy{m}(:,3),...
        'LineStyle','none',...
        'LineWidth',1,...
        'Marker','s',...
        'MarkerSize',6,...
        'Color',[60 179 113]/255,...
        'MarkerEdgeColor',[60 179 113]/255)
    
    xlabel('r (mm)')        % Radial position
    ylabel('Stresses (MPa)')% Cauchy stress
    legend('$\sigma_r$, HyperSym','$\sigma_r$, ANSYS',...
        '$\sigma_{\theta}$, HyperSym','$\sigma_{\theta}$, ANSYS',...
        '$\sigma_a$, HyperSym','$\sigma_a$, ANSYS',...
        'Location','northeast')
    
    % Save plot(s)
    GraphName = 'Stress_Distribution_Radial_Path';
    SaveFIG(h2{m},[GraphName sprintf('_%s_%s',MID{m},IntName)],subfolder,'.fig');
end
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
function [lbd,Disp] = ReadAnsysDataII(fileName,lineDisp)
% ReadAnsysDataII reads an output file generated by Ansys.

DispID = fopen(fileName,'r');

item = textscan(DispID,'%f %f %f %f',...
    'Delimiter','\n','Headerlines',lineDisp-1);
lbd = item{1};
Disp = [item{:,2} item{:,3} item{:,4}];
fclose(DispID);
end