% Cart2Cyl ----------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function Sigma = Cart2Cyl(Stress,Nodes)
% Convert a matrix in Cartesian coordinates system to cylindrical system.

Coord = Nodes(1:length(Stress),:);
theta = atan(Coord(:,2)./Coord(:,1));

Sigma = zeros(size(Stress));

for n = 1:size(Stress,2)
    % Stress vector in Cartesian coord. system
    Sv = Stress(:,n);
    
    M = [ cos(theta(n)) sin(theta(n)) 0;
         -sin(theta(n)) cos(theta(n)) 0;
                      0             0 1];
    
    % Stress matrix in Cartesian coord. system
    Sxyz = [...
        Sv(1,1) Sv(4,1) Sv(6,1);
        Sv(4,1) Sv(2,1) Sv(5,1);
        Sv(6,1) Sv(5,1) Sv(3,1)];
    
    % Stress matrix in cylindrical coord. system
    Srtz = M*Sxyz*M';
    
    % Stress vector in cylindrical coord. system
    Sigma(:,n) = [...
        Srtz(1,1); Srtz(2,2); Srtz(3,3);
        Srtz(1,2); Srtz(2,3); Srtz(1,3)];
end
end