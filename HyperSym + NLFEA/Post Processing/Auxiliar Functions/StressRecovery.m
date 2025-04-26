% Directory ---------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function [NodalStress,ElemStress] = StressRecovery(CAUCHY,XYZ,LE)
% Performs a local stress extrapolation based on the least squares theory. 
%
% Reference:
% [1] Hinton, E., & Campbell, J.S. (1974). Local and global smoothing of 
%     discontinuous finite element functions using a least squares method. 
%     International Journal for Numerical Methods in Engineering, 8(3), 
%     461-480.

TMPCAUCHY = zeros(6,size(XYZ,1));
TMP2CAUCHY = zeros(6,size(LE,1));
count = zeros(1,size(XYZ,1));

for IE = 1:size(LE,1)
    ELXY = XYZ(LE(IE,:),:);
    ElemCAUCHY = CAUCHY(:,(IE-1)*8+1:8*IE);
    
    % Extrapolation to the nodes in the full integration
    S = GaussPt2Node(ElemCAUCHY,ELXY);
    
    % Store all converged PK2 and Cauchy stresses
    TMPCAUCHY(:,LE(IE,:)) = TMPCAUCHY(:,LE(IE,:)) + S';
    
    count(:,LE(IE,:)) = count(:,LE(IE,:)) + 1;
    
    TMP2CAUCHY(:,IE) = (mean(S,1))';
end

NodalStress = TMPCAUCHY./count;   ElemStress = TMP2CAUCHY;
end
%% NODAL STRESSES
function S = GaussPt2Node(CAUCHY,ELXY)
% Integration points
XG = [-0.57735026918963D0, 0.57735026918963D0];

% Shape fnc. evaluated in the nodal pts.:
% In the shape function matrix N,
% -> each row is the value of a shape function Ni, for i = 1, 2, ..., 8 and
% -> each column corresponds to a Gauss' point.

N = zeros(8,8); line = 1;

while line <= length(N)
    for LX = 1:2
        for LY = 1:2
            for LZ = 1:2
                E1 = XG(LX); E2 = XG(LY); E3 = XG(LZ);
                
                % Determinant and shape function derivatives
                SF = SHAPEL([E1 E2 E3],ELXY);     N(line,:) = SF';
                
                line = line + 1;
            end
        end
    end
end

A = N'*N;
B = N*CAUCHY';   S = A\B;
end 