% HyperFit_Comparison_neoHookean ------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% 
% Description: This code plots the result obtained in HyperFitScript with 
% the theorethical neo-Hookean model for the test data from [1].
%
% References:
% [1] Treloar, L.R.G. (1944). Stress-strain data for vulcanized rubber 
%     under various types of deformation. Rubber Chemistry and Technology, 
%     17(4), 813-825.
% [2] Treloar, L.R.G. (1975). The Physics of Rubber Elasticity, 3rd 
%     Edition. Clarendon Press.
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
%% Compute nominal stress for the models
% Fitted model: fitted model
ua_fit = t_uniaxial(ua_exp(:,1),p);
eb_fit = t_equibiaxial(eb_exp(:,1),p);
ps_fit = t_pureshear(ps_exp(:,1),p);
% neo-Hookean: theorethical approximation from initial shear modulus
G = 0.39; % Shear modulus [2]
A10 = G/2; % From consistency conditions
ua_nH = (2*A10*(ua_exp(:,1).^3 - 1))./ua_exp(:,1).^2;
eb_nH = (2*A10*(eb_exp(:,1).^6 - 1))./eb_exp(:,1).^5;
ps_nH = (2*A10*(ps_exp(:,1).^4 - 1))./ps_exp(:,1).^3;
%% First plot
figure('Name','Compare-nH','NumberTitle','off','WindowState','maximized')
subplot(2,2,1); % title('All Tests')
hold on
plot(ua_exp(:,1),ua_exp(:,2),...
    'Color',[220 20 60]/255,'Marker','o','LineStyle','none')
plot(eb_exp(:,1),eb_exp(:,2),...
    'Color',[60 179 113]/255,'Marker','^','LineStyle','none')
plot(ps_exp(:,1),ps_exp(:,2),'bs',...
    'Color',[0 0 1],'Marker','s','LineStyle','none')
plot(ua_exp(:,1),ua_fit,...
        'Color',[220 20 60]/255,'Marker','none','LineStyle','-')
plot(eb_exp(:,1),eb_fit,...
        'Color',[60 179 113]/255,'Marker','none','LineStyle','-')
plot(ps_exp(:,1),ps_fit,...
        'Color',[0 0 1],'Marker','none','LineStyle','-')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Nominal Stress (MPa)')
leg = [...
    "Uniaxial (exp)","Equi-biaxial (exp)","Pure shear (exp)"; ...
    "Uniaxial (fit)","Equi-biaxial (fit)","Pure shear (fit)"]';
legend(leg,'Location','nw','NumColumns',3,'Orientation','horizontal')
%% Add subplot: Uniaxial test
subplot(2,2,2); % title('Simple Tension')
hold on
plot(ua_exp(:,1),ua_exp(:,2),...
    'Color',[220 20 60]/255,'Marker','o','LineStyle','none')
plot(ua_exp(:,1),ua_fit,...
    'Color',[220 20 60]/255,'Marker','none','LineStyle','-')
plot(ua_exp(:,1),ua_nH,...
        'Color',[220 20 60]/255,'Marker','none','LineStyle','-.')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Nominal Stress (MPa)')
legend(["Experiment","Ogden","neo-Hookean"],'Location','best')
%% Add subplot: Equi-biaxial test
subplot(2,2,3); % title('Equi-biaxial Tension')
hold on
plot(eb_exp(:,1),eb_exp(:,2),...
    'Color',[60 179 113]/255,'Marker','^','LineStyle','none')
plot(eb_exp(:,1),eb_fit,...
    'Color',[60 179 113]/255,'Marker','none','LineStyle','-')
plot(eb_exp(:,1),eb_nH,...
    'Color',[60 179 113]/255,'Marker','none','LineStyle','-.')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Nominal Stress (MPa)')
legend(["Experiment","Ogden","neo-Hookean"],'Location','best')
%% Add subplot: Pure shear test
subplot(2,2,4); % title('Pure Shear')
hold on
plot(ps_exp(:,1),ps_exp(:,2),...
    'Color',[0 0 1],'Marker','s','LineStyle','none')
plot(ps_exp(:,1),ps_fit,...
    'Color',[0 0 1],'Marker','none','LineStyle','-')
plot(ps_exp(:,1),ps_nH,...
        'Color',[0 0 1],'Marker','none','LineStyle','-.')
hold off; grid on; box on; axis tight
xlabel('$\lambda_1$'); ylabel('Nominal Stress (MPa)')
legend(["Experiment","Ogden","neo-Hookean"],'Location','best')