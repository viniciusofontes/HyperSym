% HYPER3D_SRI -------------------------------------------------------------
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
function [INTFORCE,GKF,SIGMA] = ...
    HYPER3D_SRI(MID,PROP,UPDATE,LTAN,XYZ,LE,DISPTD,DSFS)
% *************************************************************************
% MAIN PROGRAM FOR HYPERELASTIC MATERIAL MODELS
% COMPUTE GLOBAL STIFFNESS MATRIX RESIDUAL FORCE
%
% WARNING: This versions of HYPER3D uses the Selective Reduced Integration
% method (SRI), in which the isochoric term is computed using full
% integration (all Gauss' points) but the volumetric part uses reduced
% integration (only 1 Gauss' point)
% *************************************************************************

INTN = 0;                 % Index for history variables (each integ. pt.)
[NUMNP,NDOF] = size(XYZ); % Number of nodes and DOFs
[NE,NNE] = size(LE);
NDof = 3*NNE;             % Number of Dofs per element
NKE = NDof^2;
GKF = zeros(NE*NDof^2,1); % Global stiffness matrix
NEQ = NDOF*NUMNP;         % Number of equations (total DOFs)
INTFORCE = zeros(NEQ,1);  % Int. force vectors
NDOFE = NDOF*8;           % Number of DOFs per element
SIGMA = zeros(6,8*NE);
indexK = 0;

% LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
for IE = 1:NE
    EKF = zeros(NDOFE,NDOFE);

    % Nodal coordinates and incremental displacements
    ELXY = XYZ(LE(IE,:),:);

    % Local to global mapping
    IDOF = zeros(1,24);

    for I = 1:8
        II = (I-1)*NDOF+1;
        IDOF(II:II+2) = (LE(IE,I)-1)*NDOF+1:(LE(IE,I)-1)*NDOF+3;
    end

    DSP = DISPTD(IDOF);
    DSP = reshape(DSP,NDOF,8);

    % LOOP OVER INTEGRATION POINTS (only for the isochoric part)
    % Integration points and weights for the isochoric part
    % (full integration)    
    WGT = [1.00000000000000D0, 1.00000000000000D0];

    for LX = 1:2
        for LY = 1:2
            for LZ = 1:2
                INTN = INTN + 1;

                % Determinant and shape function derivatives
                INDEX = LX*4 + LY*2 + LZ - 6;
                DSF  = DSFS(:,:,INDEX);
                GJ   = DSF*ELXY;
                DET  = det(GJ);
                SHPD = GJ\DSF;
                FAC  = WGT(LX)*WGT(LY)*WGT(LZ)*DET;

                % Deformation gradient
                F = DSP*SHPD' + eye(3);

                % Right Cauchy-Green tensor
                C = F'*F;

                % Computer stress and tangent stiffness
                STRESS = zeros(6,1);
                STRESS(:) = feval([MID '_iso'],C,PROP);

                if UPDATE % Update Cauchy stress
                    SIGMA(:,INTN) = CAUCHY(F,STRESS);
                    continue;
                end

                % Residual force and tangent stiffness matrix
                BN = zeros(6,24);
                BG = zeros(9,24);

                for I = 1:8
                    COL = (I-1)*3+1:(I-1)*3 + 3;

                    BN(:,COL) = [...
                        SHPD(1,I)*F(1,1)                    SHPD(1,I)*F(2,1)                    SHPD(1,I)*F(3,1);
                        SHPD(2,I)*F(1,2)                    SHPD(2,I)*F(2,2)                    SHPD(2,I)*F(3,2);
                        SHPD(3,I)*F(1,3)                    SHPD(3,I)*F(2,3)                    SHPD(3,I)*F(3,3);
                        SHPD(1,I)*F(1,2) + SHPD(2,I)*F(1,1) SHPD(1,I)*F(2,2) + SHPD(2,I)*F(2,1) SHPD(1,I)*F(3,2) + SHPD(2,I)*F(3,1);
                        SHPD(2,I)*F(1,3) + SHPD(3,I)*F(1,2) SHPD(2,I)*F(2,3) + SHPD(3,I)*F(2,2) SHPD(2,I)*F(3,3) + SHPD(3,I)*F(3,2);
                        SHPD(1,I)*F(1,3) + SHPD(3,I)*F(1,1) SHPD(1,I)*F(2,3) + SHPD(3,I)*F(2,1) SHPD(1,I)*F(3,3) + SHPD(3,I)*F(3,1)];

                    BG(:,COL) = [...
                        SHPD(1,I)         0         0;
                        SHPD(2,I)         0         0;
                        SHPD(3,I)         0         0;
                                0 SHPD(1,I)         0;
                                0 SHPD(2,I)         0;
                                0 SHPD(3,I)         0;
                                0         0 SHPD(1,I);
                                0         0 SHPD(2,I);
                                0         0 SHPD(3,I)];
                end

                % Residual and internal forces
                AUX = FAC*BN'*STRESS;
                INTFORCE(IDOF) = INTFORCE(IDOF) + AUX;

                % Tangent stiffness
                if LTAN
                    DTAN = zeros(36,1);
                    [~,DTAN(:)] = feval([MID '_iso'],C,PROP);
                    DTAN = reshape(DTAN,6,6);

                    SIG = [...
                        STRESS(1) STRESS(4) STRESS(6);
                        STRESS(4) STRESS(2) STRESS(5);
                        STRESS(6) STRESS(5) STRESS(3)];

                    SHEAD = zeros(9);
                    SHEAD(1:3,1:3) = SIG;
                    SHEAD(4:6,4:6) = SIG;
                    SHEAD(7:9,7:9) = SIG;

                    EKF = EKF + FAC*(BN'*DTAN*BN + BG'*SHEAD*BG);
                end
            end
        end
    end

    % LOOP OVER 1 INTEGRATION POINT (only for the volumetric part)
    % Integration points and weights for the isochoric part
    % (reduced integration)
    WGT = 2.00000000000000D0;

    % Determinant and shape function derivatives
    DSF  = DSFS(:,:,9);
    GJ   = DSF*ELXY;
    DET  = det(GJ);
    SHPD = GJ\DSF;
    FAC  = WGT(1)*WGT(1)*WGT(1)*DET;

    % Deformation gradient
    F = DSP*SHPD' + eye(3);

    % Right Cauchy-Green strain vector
    C = F'*F;

    % Computer stress and tangent stiffness
    STRESS = zeros(6,1);
    STRESS(:) = feval([MID '_vol'],C,PROP);

    % The STRESS onsidered in here is the volumetric part
    % (second column) only

    if UPDATE % Update Cauchy stresses
        SIGMA(:,INTN-7:INTN) = SIGMA(:,INTN-7:INTN) + CAUCHY(F,STRESS);
        continue;
    end

    % Add residual force and tangent stiffness matrix
    BN = zeros(6,24);
    BG = zeros(9,24);

    for I = 1:8
        COL = (I-1)*3+1:(I-1)*3 + 3;

        BN(:,COL) = [...
            SHPD(1,I)*F(1,1)                    SHPD(1,I)*F(2,1)                    SHPD(1,I)*F(3,1);
            SHPD(2,I)*F(1,2)                    SHPD(2,I)*F(2,2)                    SHPD(2,I)*F(3,2);
            SHPD(3,I)*F(1,3)                    SHPD(3,I)*F(2,3)                    SHPD(3,I)*F(3,3);
            SHPD(1,I)*F(1,2) + SHPD(2,I)*F(1,1) SHPD(1,I)*F(2,2) + SHPD(2,I)*F(2,1) SHPD(1,I)*F(3,2) + SHPD(2,I)*F(3,1);
            SHPD(2,I)*F(1,3) + SHPD(3,I)*F(1,2) SHPD(2,I)*F(2,3) + SHPD(3,I)*F(2,2) SHPD(2,I)*F(3,3) + SHPD(3,I)*F(3,2);
            SHPD(1,I)*F(1,3) + SHPD(3,I)*F(1,1) SHPD(1,I)*F(2,3) + SHPD(3,I)*F(2,1) SHPD(1,I)*F(3,3) + SHPD(3,I)*F(3,1)];

        BG(:,COL) = [...
            SHPD(1,I)         0         0;
            SHPD(2,I)         0         0;
            SHPD(3,I)         0         0;
                    0 SHPD(1,I)         0;
                    0 SHPD(2,I)         0;
                    0 SHPD(3,I)         0;
                    0         0 SHPD(1,I);
                    0         0 SHPD(2,I);
                    0         0 SHPD(3,I)];
    end

    % Residual and internal forces
    AUX = FAC*BN'*STRESS;
    INTFORCE(IDOF) = INTFORCE(IDOF) + AUX;

    % Tangent stiffness
    if LTAN
        DTAN = zeros(36,1);
        [~,DTAN(:)] = feval([MID '_vol'],C,PROP);
        DTAN = reshape(DTAN,6,6);

        SIG = [...
            STRESS(1) STRESS(4) STRESS(6);
            STRESS(4) STRESS(2) STRESS(5);
            STRESS(6) STRESS(5) STRESS(3)];

        SHEAD = zeros(9);
        SHEAD(1:3,1:3) = SIG;
        SHEAD(4:6,4:6) = SIG;
        SHEAD(7:9,7:9) = SIG;

        EKF = EKF + FAC*(BN'*DTAN*BN + BG'*SHEAD*BG);
    end
    
    if LTAN
        GKF(indexK+1:indexK+NKE) = EKF(:);
        indexK = indexK + NKE;
    end
end
end
%% CAUCHY STRESS
function STRESS = CAUCHY(F,S)
%***********************************************************************
% CONVERT 2ND PK STRESS INTO CAUCHY STRESS
%***********************************************************************

PK = [...
    S(1) S(4) S(6);
    S(4) S(2) S(5);
    S(6) S(5) S(3)];

ST = F*PK*F'/det(F);

STRESS = [ST(1,1); ST(2,2); ST(3,3); ST(1,2); ST(2,3); ST(1,3)];
end