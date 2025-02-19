% NLFEA -------------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
%
% Reference:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element Analysis.
%     Springer Science & Business Media, pp. 214-216.
% -------------------------------------------------------------------------
function NLFEA(ITRA,TOL,ATOL,NTOL,TIMS,MID,PROP,EXTFORCE,SDISPT,XYZ,LE,DIR,SRI,PRESS)
% *************************************************************************
% MAIN PROGRAM FOR HYPERELASTIC ANALYSIS
% *************************************************************************

fprintf('\t *** Running nonlinear analysis *** \t\n')

[NUMNP,NDOF] = size(XYZ); 	          % Number of nodes and DOFs
NEQ = NDOF*NUMNP;                     % Number of equations (total DOFs)
SOLN = zeros(NEQ,1);                  % Solution at every iteration
[DISPTD,DISPDD] = deal(zeros(NEQ,1)); % Nodal displacement & increment
ITGZONE(XYZ,LE);			          % Check element connectivity

DSFS = TAB_SHAPEL(SRI);               % Deriv. of shape fnc. (local coords.) 
[I,J] = PREP_SPARSE(LE);

% Load increments [Start End Increment InitialLoad FinalLoad]
NLOAD = size(TIMS,2);
ILOAD = 1;						      % First load increment
TIMEI = TIMS(1,ILOAD);			      % Starting time
TIMEF = TIMS(2,ILOAD);			      % Ending time
DELTA = TIMS(3,ILOAD);			      % Time increment
CUR1 = TIMS(4,ILOAD);			      % Starting load factor
CUR2 = TIMS(5,ILOAD);			      % Ending load factor
DELTA0 = DELTA;					      % Saved time increment
TIME = TIMEI;					      % Starting time
TDELTA = TIMEF - TIMEI;			      % Time interval for load step
ITOL = 1;						      % Bisection level
TARY = zeros(NTOL,1);		          % Time stamps for bisections

% Load increment loop
ISTEP = -1; POS = 1; FLAG10 = 1;

while FLAG10 == 1 				      % Solution has been converged
    FLAG10 = 0; FLAG11 = 1; FLAG20 = 1;

    CDISP = DISPTD; 				  % Store converged displacement

    if ITOL == 1 				      % No bisection
        DELTA = DELTA0;
        TARY(ITOL) = TIME + DELTA;
    else							       % Recover previous bisection
        ITOL = ITOL - 1;			       % Reduce the bisection level
        DELTA = TARY(ITOL) - TARY(ITOL+1); % New time increment
        TARY(ITOL+1) = 0;			       % Empty conv. bisection level
        ISTEP = ISTEP - 1;			       % Decrease load increment
    end

    TIME0 = TIME;					       % Save the current time

    % Update stresses and history variables
    UPDATE = true; LTAN = false;

    if ~SRI    % Full integration
        [~,~,SIGMA] = HYPER3D(MID,PROP,UPDATE,LTAN,XYZ,LE,DISPTD,DSFS);
    elseif SRI % Nearly incompressible models with SRI
        [~,~,SIGMA] = HYPER3D_SRI(MID,PROP,UPDATE,LTAN,XYZ,LE,DISPTD,DSFS);
    else
        error('Wrong material ID.');
    end

    if ISTEP >= 0
        PROUT(DIR,POS,FACTOR,DISPTD,INTFORCE,SIGMA)
        POS = POS + 1;
    end

    TIME = TIME + DELTA;		              % Increase time
    ISTEP = ISTEP + 1;

    % Check time and control bisection
    while FLAG11 == 1 			              % Bisection loop start
        FLAG11 = 0;

        if TIME - TIMEF  > 1E-10 	          % Time passed the end time
            if TIMEI + DELTA - TIME  > 1E-10  % One more at the end time
                DELTA = TIMEF + DELTA - TIME; % Time increment to the end
                DELTA0 = DELTA;			      % Saved time increment
                TIME = TIMEF;			      % Current time is the end
            else
                ILOAD = ILOAD + 1;		      % Progress to next load step
                if ILOAD > NLOAD 		      % Finished final load step
                    FLAG10 = 0;			      % Stop the program
                    break;
                else					      % Next load step
                    TIME = TIME - DELTA;
                    DELTA = TIMS(3,ILOAD);
                    DELTA0 = DELTA;
                    TIME = TIME + DELTA;
                    TIMEI = TIMS(1,ILOAD);
                    TIMEF = TIMS(2,ILOAD);
                    TDELTA = TIMEF - TIMEI;
                    CUR1 = TIMS(4,ILOAD);
                    CUR2 = TIMS(5,ILOAD);
                end
            end
        end

        % Load factor and prescribed displacements
        FACTOR = CUR1 + (TIME - TIMEI)/TDELTA*(CUR2 - CUR1);
        SDISP = DELTA*SDISPT(:,3)/TDELTA*(CUR2 - CUR1);

        % Start convergence iteration
        ITER = 0;

        while FLAG20 == 1
            FLAG20 = 0;
            ITER = ITER + 1;

            % Assemble K and F
            UPDATE = false; LTAN = true;

            if ~SRI    % Nearly incompressible models with full integration
                [INTFORCE,GKF] = HYPER3D(MID,PROP,UPDATE,LTAN,XYZ,LE,DISPTD,DSFS);
            elseif SRI % Nearly incompressible models with SRI
                [INTFORCE,GKF] = HYPER3D_SRI(MID,PROP,UPDATE,LTAN,XYZ,LE,DISPTD,DSFS);
            end

            GKF = sparse(I,J,GKF);

            % Increase external force
            FORCE = -INTFORCE;

            if size(EXTFORCE,1) > 0
                if exist('PRESS','var') == 1
                    EXTFORCE = PressureLoad(XYZ,LE,PRESS.EREF,PRESS.Pi,DISPTD);
                end

                LOC = NDOF*(EXTFORCE(:,1)-1) + EXTFORCE(:,2);
                FORCE(LOC) = FORCE(LOC) + FACTOR*EXTFORCE(:,3);
            end

            % Prescribed displacement BC
            NDISP = size(SDISPT,1);

            if NDISP ~= 0
                FIXEDDOF = NDOF*(SDISPT(:,1)-1) + SDISPT(:,2);
                GKF(FIXEDDOF,:) = zeros(NDISP,NEQ);
                GKF(FIXEDDOF,FIXEDDOF) = PROP(1)*eye(NDISP);

                FORCE(FIXEDDOF) = 0;

                if ITER == 1, FORCE(FIXEDDOF) = PROP(1)*SDISP; end
            end

            % Check convergence
            if ITER > 1
                FIXEDDOF = NDOF*(SDISPT(:,1)-1) + SDISPT(:,2);
                ALLDOF = 1:NEQ;
                FREEDOF = setdiff(ALLDOF,FIXEDDOF);

                if size(SDISPT,1) < NEQ
                    RESN = max(abs(FORCE(FREEDOF)));
                else
                    RESN = norm(SOLN);
                end

                OUTPUT(1,ITER-1,RESN,TIME,DELTA)

                if RESN < TOL, FLAG10 = 1; break; end

                if  RESN > ATOL || ITER >= ITRA       % Start bisection
                    ITOL = ITOL + 1;

                    if ITOL <= NTOL
                        DELTA = 0.5*DELTA;
                        TIME = TIME0 + DELTA;
                        TARY(ITOL) = TIME;
                        DISPTD = CDISP;
                        fprintf(1,['\n   Not converged!\n'...
                            '   Bisecting load increment %2d\n'],ITOL-1)
                    else
                        fprintf(1,'\t\t *** Max No. of bisection ***\n')
                        return;
                    end

                    FLAG11 = 1; FLAG20 = 1; break;
                end
            end

            % Solve the system equation
            if FLAG11 == 0
                SOLN = GKF\FORCE;
                DISPDD = DISPDD + SOLN;  DISPTD = DISPTD + SOLN;
                FLAG20 = 1;
            else
                FLAG20 = 0;
            end

            if FLAG10 == 1, break; end
        end 							% 20 Convergence iteration
    end 								% 11 Bisection
end 								    % 10 Load increment
end
%% CHECKING CONVERGENCE
function VOLUME = ITGZONE(XYZ,LE)
%*************************************************************************
% Check element connectivity and calculate volume
%*************************************************************************

EPS = 1E-7;
NE = size(LE,1);
VOLUME = 0;

for I = 1:NE
    ELXY = XYZ(LE(I,:),:);
    [~,~,DET] = SHAPEL([0 0 0],ELXY);
    DVOL = 8*DET;

    if DVOL < EPS, error('Negative Jacobian'); end

    VOLUME = VOLUME + DVOL;
end
end
%% OUTPUT CONVERGENCE HISTORY
function OUTPUT(FLG,ITER,RESN,TIME,DELTA)
%*************************************************************************
% Print convergence iteration history
%*************************************************************************

if FLG == 1
    if ITER > 1
        fprintf(1,'%27d  %10.5e \n',ITER,full(RESN))
    else
        fprintf(1,'\n \t  Time  Time step  Iter  Residual \n')
        fprintf(1,'%10.5f %10.3e %5d  %10.5e \n',TIME,DELTA,ITER,full(RESN))
    end
end
end
%% PRINTING VARIABLES IN '.MAT' FILE
function PROUT(DIR,t,FACTOR,DISPTD,INTFORCE,SIGMA)
%*************************************************************************
% Print converged displacements, load factor and stresses in the current
% step
%*************************************************************************

if isfile(DIR)
    load(DIR,'FAC','DISP','Fi','CAUCHY');
end
    
% Store converged nodal displacements
DISP(:,:,t) = (reshape(DISPTD,3,[]))';

% Store current load factor
FAC(t) = FACTOR;

% Store internal force
Fi(:,:,t) = reshape(full(INTFORCE),3,[])';

% Store cauchy stress at Gauss points
CAUCHY(:,:,t) = SIGMA;
pause(.1)
save(DIR,'FAC','DISP','Fi','CAUCHY');
end