% HyperFit_Relative_Error -------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
%
% Description: This code presents plots with the relative errors for the
% result obtained in HyperFitScript
% -------------------------------------------------------------------------
%% Set defaults
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureColor',[1 1 1])
set(0,'defaultAxesFontSize',18)
set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultLineLineWidth',2);
%% First plot 
figure('Name','Relative Error','NumberTitle','off','WindowState','maximized')
subplot(2,2,1); % title('All Tests')
hold on
plot(Error(1).rel(:,1),100*Error(1).rel(:,2),...
    'Color',[220 20 60]/255,'Marker','x','LineStyle','-')
plot(Error(2).rel(:,1),100*Error(2).rel(:,2),...
    'Color',[60 179 113]/255,'Marker','x','LineStyle','-')
plot(Error(3).rel(:,1),100*Error(3).rel(:,2),...
    'Color',[0 0 1],'Marker','x','LineStyle','-')
plot(Error(1).rel(:,1),5*ones(size(Error(1).rel,1)),...
    'Color',[0 0 0],'Marker','none','LineStyle','--')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Relative Error (\%)')
leg = ["Uniaxial","Equi-biaxial","Pure shear","5\% error"]';
legend(leg,'Location','nw','NumColumns',4,'Orientation','horizontal')
%% Add subplot: Uniaxial test
subplot(2,2,2); % title('Simple Tension')
hold on
plot(Error(1).rel(:,1),100*Error(1).rel(:,2),...
    'Color',[220 20 60]/255,'Marker','x','LineStyle','-')
plot(Error(1).rel(:,1),5*ones(size(Error(1).rel,1)),...
    'Color',[0 0 0],'Marker','none','LineStyle','--')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Relative Error (\%)')
legend(["Fitted model","5\% error"],'Location','ne')
%% Add subplot: Equi-biaxial test
subplot(2,2,3); % title('Equi-biaxial Tension')
hold on
plot(Error(2).rel(:,1),100*Error(2).rel(:,2),...
    'Color',[60 179 113]/255,'Marker','x','LineStyle','-')
plot(Error(2).rel(:,1),5*ones(size(Error(2).rel,1)),...
    'Color',[0 0 0],'Marker','none','LineStyle','--')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Relative Error (\%)')
legend(["Fitted model","5\% error"],'Location','ne')
%% Add subplot: Pure shear test
subplot(2,2,4); % title('Pure Shear')
hold on
plot(Error(3).rel(:,1),100*Error(3).rel(:,2),...
    'Color',[0 0 1],'Marker','x','LineStyle','-')
plot(Error(3).rel(:,1),5*ones(size(Error(3).rel,1)),...
    'Color',[0 0 0],'Marker','none','LineStyle','--')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Relative Error (\%)')
legend(["Fitted model","5\% error"],'Location','ne')