% TAB_SHAPEL---------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
%
% Reference:
% [1] Kim, N.-H. (2014). Introduction to Nonlinear Finite Element Analysis.
%     Springer Science & Business Media, pp. 214-216.
% -------------------------------------------------------------------------
function DSFS = TAB_SHAPEL(SRI)
% *************************************************************************
% COMPUTE AND STORE SHAPE FUNCTIONS (DERIVATIVES AND JACOBIANS)
% *************************************************************************

% Integration points and weights
XG = [-0.57735026918963D0, 0.57735026918963D0];

% Number of integration points (2x2x2=8, +1 at centroid if SRI enabled)
if SRI, NINT = 9; else, NINT = 8; end

% Derivatives of shape functions (in the local coordinates XI)
DSFS = zeros(3,8,NINT);

XNODE = [...
    -1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
QUAR = 0.125;

% LOOP OVER ELEMENTS, THIS IS MAIN LOOP TO COMPUTE K AND F
for LX = 1:2
    for LY = 1:2
        for LZ = 1:2
            DSF = zeros(3,8);
            XI = [XG(LX) XG(LY) XG(LZ)];

            for I = 1:8
                XP = XNODE(1,I);
                YP = XNODE(2,I);
                ZP = XNODE(3,I);

                XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];

                DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
                DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
                DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
            end

            INDEX = LX*4 + LY*2 + LZ - 6;
            DSFS(:,:,INDEX) = DSF;
        end
    end
end

if SRI
    DSF = zeros(3,8);
    XI = [0 0 0];

    for I = 1:8
        XP = XNODE(1,I);
        YP = XNODE(2,I);
        ZP = XNODE(3,I);

        XI0 = [1+XI(1)*XP 1+XI(2)*YP 1+XI(3)*ZP];

        DSF(1,I) = QUAR*XP*XI0(2)*XI0(3);
        DSF(2,I) = QUAR*YP*XI0(1)*XI0(3);
        DSF(3,I) = QUAR*ZP*XI0(1)*XI0(2);
    end

    DSFS(:,:,9) = DSF;
end
end