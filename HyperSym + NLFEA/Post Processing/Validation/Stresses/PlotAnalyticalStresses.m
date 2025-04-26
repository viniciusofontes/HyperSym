% PlotAnalyticalStresses -------------------------------------------------- 
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function [h1,h2] = PlotAnalyticalStresses(E0,nu0,MID)
% Symbolic variables
syms lambda0 mu0 gamma0 real

% Compute kinematic variables
I = sym(eye(3));        % Identity tensor
F = I; F(1,1) = gamma0; % Deformation gradient
C = F'*F;               % Right Cauchy-Green strain
Cinv = inv(C);          % Inverse of C
E = (C - I)/2;          % Green-Lagrange strain tensor
J = det(F);             % Jacobian of deformation
%% Analytical expressions for PK2 (Klarbring & Strömberg, 2013)
PK2{1} = [lambda0*trace(E)*I 2*mu0*E];
PK2{2} = [lambda0*log(J)./C 2*mu0*E];
PK2{3} = [lambda0*(J - 1)./C 2*mu0*E];
PK2{4} = [lambda0*(J - 1)*J./C 2*mu0*E];
PK2{5} = [lambda0*log(J)./C  mu0*(I - Cinv)];
PK2{6} = [lambda0*(J - 1)./C  mu0*(I - Cinv)];
PK2{7} = [lambda0*(J - 1)*J./C  mu0*(I - Cinv)];

% Iterate through models and push-forward PK2
sigma11 = sym(zeros(7,2));

fprintf('Analytical Stresses =========================\n\n')

for m = 1:7
    cauchyL = 1/J*F*PK2{m}(:,1:3)*F';  cauchyR = 1/J*F*PK2{m}(:,4:6)*F';
    sigma11(m,:) = simplify([cauchyL(1,1) cauchyR(1,1)]);
    
    fprintf('sigma_%d(1,1) = \n',m); pretty((sigma11(m,:)))
end

Lambda = E0*nu0/((1 + nu0)*(1 - 2*nu0));
Mu = E0/(2*(1 + nu0));
Gamma = linspace(0.2,2)';

S = zeros(length(Gamma),7);

for m = 1:7
    S(:,m) = eval(subs(sum(sigma11(m,:)),{gamma0,lambda0,mu0},...
        {Gamma,Lambda,Mu}));
end

Color = {[1 0.8 0]; [0 0 1]; [0 0 1]; [0 0 1]; [1 0 0]; [1 0 0]; [1 0 0]};
LineStyle = {'-','-','-.','--','-','-.','--'};
%% Plot
% SVH-based material models -----------------------------------------------
MinS = min(min(S(:,1:4)));
MaxS = max(max(S(:,1:4)));

MinRange = 0.9;
MaxRange = 2 - MinRange;

h1 = figure; hold all; axis tight; box on

patch('Faces',1:4,...
    'Vertices',[MinRange MinS; MaxRange MinS; MaxRange MaxS; MinRange MaxS],...
    'FaceAlpha',0.1,...
    'EdgeColor','none',...
    'HandleVisibility','off')

for i = 1:4
    plot(Gamma,S(:,i),...
        'Color',Color{i},...
        'LineStyle',LineStyle{i},...
        'LineWidth',1,...
        'DisplayName',['Analytical, ' MID{i}])
end

xticks(sort([0.2:0.2:2 MinRange MaxRange]))
xlabel('$\gamma$'); ylabel('$\sigma_{11}$')

legend('-DynamicLegend','Location','southeast')

% nH-based material models ------------------------------------------------
MinS = min(min(S(:,5:7)));
MaxS = max(max(S(:,5:7)));

% MinRange = 0.9;
% MaxRange = 2 - MinRange;

h2 = figure; hold all; axis tight; box on

patch('Faces',1:4,...
    'Vertices',[MinRange MinS; MaxRange MinS; MaxRange MaxS; MinRange MaxS],...
    'FaceAlpha',0.1,...
    'EdgeColor','none',...
    'HandleVisibility','off')

for j = 5:7
    plot(Gamma,S(:,j),...
        'Color',Color{j},...
        'LineStyle',LineStyle{j},...
        'LineWidth',1,...
        'DisplayName',['Analytical, ' MID{j}])
end

xticks(sort([0.2:0.2:2 MinRange MaxRange]))
xlabel('$\gamma$'); ylabel('$\sigma_{11}$')

legend('-DynamicLegend','Location','southeast')
end