% HyperFit ----------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function [p_opt,Error,res_opt] = HyperFit(W,PROP,Test,opt)
% Fit parameters from experimental data and a given material model.
[NTest,NPROP] = deal(length(Test),length(PROP)); % # of tests and props.
% Apply inconpressibility condition (J = lambda1*lambda2*lambda3 = 1)
syms lambda1 lambda2 lambda3
W = subs(W,lambda3,1/(lambda1*lambda2));
%% Derive analytical expressions for the nominal stress in each test
% Initialize nominal stress (t) and its gradient w.r.t. properties
Dir = [cd '/Functions/'];  % Directory to save the functions
if ~exist(Dir,'dir'), mkdir(Dir); end
for n = 1:NTest
    t = feval(Test(n).name,W); % Nominal stress for a given test
    dtdPROP = simplify(gradient(t,PROP));
    % Generate MATLAB functions
    Controls = struct('File',[Dir 't_' Test(n).name],'Optimize',true); 
    matlabFunction(t,dtdPROP,'Vars',{lambda1,PROP},Controls);
end
%% Minimize square error to find material parameters
npoints = zeros(length(Test),1);
[xdata,ydata] = deal(zeros(0,1));
for n = 1:NTest
    npoints(n) = size(Test(n).data,1);
    xdata = cat(1,xdata,Test(n).data(:,1));
    ydata = cat(1,ydata,Test(n).data(:,2));
end
% Create obj. fnc.: x = mat. parameters; xdata = principal stretch values
fun = @(x,xdata) ObjFnc(x,xdata,Test,npoints);
options = optimoptions('lsqcurvefit',...
    'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',1e-12,...
    'OptimalityTolerance',1e-12,...
    'FunctionTolerance',1e-12,...
    'MaxIterations',1e3,...
    'MaxFunctionEvaluations',1e6,...
    'Display','none');
% Set parameters to default or value in opt (if provided)
if nargin < 4, opt = struct(); end % Initialize opt if not provided
N  =  setOpt(opt,'N',10);  % Number of optimization runs
lb =  setOpt(opt,'pL',[]); % Lower boundary of p
ub =  setOpt(opt,'pU',[]); % Upper boundary of p
lb0 = setOpt(opt,'pL',-1); % Lower boundary of initial guess
ub0 = setOpt(opt,'pU',+1); % Upper boundary of initial guess
x0_range = ub0 - lb0;      % Range of initial guess
[p_opt,x0] =  deal(setOpt(opt,'p0',[])); % First initial guess
% Make boundaries a uniform vector if it is a scalar
if length(lb) == 1, lb = lb*ones(size(x0)); end
if length(ub) == 1, ub = ub*ones(size(x0)); end
% Run optimization for N initial guesses
for n = 1:N 
    fprintf('Curve fit #%d:\t',n)
    if n ~= 1 || isempty(x0) % Initial guess unless provided and n == 1
        x0 = x0_range*rand(1,NPROP) + lb0;
    end
    [p,resnorm] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
    if n == 1 || resnorm < res_opt
        fprintf('New optimum found!\n')
        [n_opt,p_opt,res_opt] = deal(n,p,resnorm);
    else
        fprintf('\n')
    end
end
if N > 0 % Further optimize the best result with tighter tolerances    
    [options.StepTolerance,options.OptimalityTolerance] = deal(1e-15);
    options.MaxIterations = 1e4; % Increase optimization iterations
    fprintf('\nFurther optimizing the solution #%d\n',n_opt);
    [p_opt,res_opt] = lsqcurvefit(fun,p_opt,xdata,ydata,lb,ub,options);
end
%% Check consistency condition (c.c.)
cc = setOpt(opt,'cc',false); % Derive and check c.c.?
if cc % Estimate shear modulus from c.c.    
    dWdl1 = diff(W,lambda1);
    dWdl1dl1 = diff(dWdl1,lambda1);
    dWdl2dl1 = diff(dWdl1,lambda2);
    lhs = dWdl1 + dWdl1dl1 - dWdl2dl1; % Left-hand side of c.c.   
    lhs = simplify(subs(lhs,[lambda1 lambda2],[1 1]));
    shear_modulus_cc = vpa(subs(lhs,PROP,p_opt))/2;
    fprintf('\nConsistency condition (c.c.) for the current model\n\n')
    syms mu1, pretty(simplify(lhs == 2*mu1))
    fprintf('Shear modulus from c.c. and optimum parameters: %.3f\n\n',shear_modulus_cc)    
end
%% Compute error and relative error for the optimal set of properties
if nargout > 1
    Stress = ObjFnc(p_opt,xdata,Test,npoints);
    Error = struct([]);    
    [b,res_opt] = deal(0);
    for n = 1:NTest
        a = b + 1;              % Index of first data point in the curr. test
        b = a + npoints(n) - 1; % Index of last data point of curr. test
        Error(n).abs  = [xdata(a:b) abs(Stress(a:b) - ydata(a:b))];
        Error(n).quad = sum((Error(n).abs(:,2)).^2);
        Error(n).rel  = [xdata(a:b) Error(n).abs(:,2)./(max(0.5,abs(ydata(a:b))))];    
        res_opt = res_opt + Error(n).quad;
    end
end
end
%% Uniaxial test
function t1 = uniaxial(W)
% Kinematic constraints for each uniaxial test.
syms lambda1 lambda2
W = subs(W,lambda2,1/sqrt(lambda1)); % Apply kinematic constraint
t1 = simplify(diff(W,lambda1));       % t = dW/dlambda1
end
%% Equi-biaxial test
function t1 = equibiaxial(W)
% Kinematic constraints for each equi-biaxial test.
syms lambda1 lambda2
W = subs(W,lambda2,lambda1);      % Apply kinematic constraint
t1 = simplify(diff(W,lambda1))/2; % t = dW/dlambda1/2 (lambda1 = lambda2)
end
%% (Pure) shear test
function t1 = pureshear(W)
% Kinematic constraints for each pure shear test.
syms lambda1 lambda2
W = subs(W,lambda2,1);          % Apply kinematic constraint
t1 = simplify(diff(W,lambda1)); % t = dW/dlambda1
end
%% Set optimization parameters
function value = setOpt(opt,field,default)
% Retrieve optimization parameters from opt.
if isfield(opt,field) && ~isempty(opt.(field))
    value = opt.(field);
else
    value = default;
end
end
%% Objective function
function [f,g] = ObjFnc(p,xdata,Test,npoints)
% Substitutive the current value of the properties in the expressions for
% the objective function.
[NTest,NPROP,Ndata] = deal(length(Test),length(p),length(xdata));
f = zeros(size(xdata)); g = zeros(Ndata,NPROP); b = 0;
for n = 1:NTest
    a = b + 1;              % Index of first data point in the curr. test
    b = a + npoints(n) - 1; % Index of last data point of curr. test
    lambda1 = xdata(a:b);   % Principal stretches
    % Nominal Stress (from analytical model)
    [f(a:b),gn] = feval(['t_' Test(n).name],lambda1,p);
    g(a:b,:) = reshape(gn,[npoints(n) NPROP]);
end
end