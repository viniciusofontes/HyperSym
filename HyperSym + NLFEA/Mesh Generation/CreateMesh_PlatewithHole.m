% CreateMesh_PlatewithHole ------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2025). 
%          HyperSym: an educational MATLAB code for hyperelasticity
%          Computer Applications in Engineering Education
%          DOI: 10.1002/cae.70037
% -------------------------------------------------------------------------
function [Nodes,Elem] = CreateMesh_PlatewithHole(Dim,BdBox,NElem)
% Generate mesh. It starts from a 2D model which is extruded.
% Then the element connectivities are performed.
%
% OBS.: 
% The number of elements, NElem, corresponds to half of geometry 
% considering the symmetry of line y = x

NElemX = NElem(1);    % number of elements along X-axis
NElemY = NElem(2);    % number of elements along Y-axis
NElemZ = NElem(3);    % number of elements along Z-axis
          
DeltaX = (BdBox(2) - BdBox(1))/NElemX;
DeltaY = (BdBox(4) - BdBox(3))/NElemY;
DeltaZ = (BdBox(6) - BdBox(5))/NElemZ;

alpha = 90;

fprintf('Generating mesh...\n\n')

% 2D Mesh -----------------------------------------------------------------
% Node Assembling
x1 = (BdBox(1):DeltaX:BdBox(2))';
x1 = repmat(x1,1,NElemY + 1);
y1 = zeros(NElemY + 1,NElemY + 1);
i = 1; theta = 0;

while theta < 45
    i = i+1;
    theta = theta + alpha/(NElemY+NElemX);
    
    % 1st coordinate follows the inner radious. 
    % Other are fixed in x-axis or y-axis.
    x1(1,i) = Dim.r*cosd(theta);
    y1(1,i) = Dim.r*sind(theta);
    
    y1(2:end,i) = tand(theta)*x1(2:end,i);  % line y = m*x, m = atan(y/x)
end

y2 = (BdBox(3):DeltaY:BdBox(4))';
y2 = repmat(y2,1,NElemX + 1);
x2 = zeros(NElemX + 1,NElemX + 1);
j = 0;

while theta > 0
    j = j+1;
    theta = theta - alpha/(NElemY+NElemX);
    
    % 1st coordinate follows the inner radious. 
    % Other are fixed in x-axis or y-axis.
    x2(1,j) = Dim.r*sind(theta);
    y2(1,j) = Dim.r*cosd(theta);

    x2(2:end,j) = tand(theta)*y2(2:end,j);
end

% Eliminating duplicate coordinates and node assembling
x = [x1(:); x2(:)];     y = [y1(:); y2(:)];
NodesInPlane = unique([x y],'rows','stable');

% Element Assembling
ElemInPlane = cell(2*NElemX*NElemY,1);

for elx = 1:NElemX
    for ely = 1:2*NElemY
        ElemInPlane{NElemX*(ely - 1) + elx} = [...
            (NElemX + 1)*(ely - 1) + elx ...
            (NElemX + 1)*(ely - 1) + elx + 1 ...
            (NElemX + 1)*ely + elx + 1 ...
            (NElemX + 1)*ely + elx];
    end
end

% 3D Mesh -----------------------------------------------------------------
% Node Assembling
Nodes = zeros((NElemX + 1)*(2*NElemY + 1)*(NElemZ + 1),3);
Nodes(:,1:2) = repmat(NodesInPlane,NElemZ+1,1); % extrusion in z-direction

row = 0;     t = 0;     Ne = length(NodesInPlane);

while row <= NElemZ     % fill the extruded coordinates
    Nodes(row*Ne + 1:(row+1)*Ne,3) = t;
    t = t + DeltaZ;     row = row + 1;
end

% Element Assembling
NNodesInPlane = length(NodesInPlane);
Elem = cell(2*NElemX*NElemY*NElemZ,1);
Ne = length(ElemInPlane);
ElemI = ElemInPlane; ElemJ = ElemI;

for elz = 1:NElemZ
    aux = 1;
    
    for z = 1+(elz-1)*Ne:elz*Ne
        ElemJ{aux} = ElemInPlane{aux} + elz*NNodesInPlane;
        Elem{z} = [ElemI{aux} ElemJ{aux}];
        
        aux = aux + 1;
    end
    
    ElemI = ElemJ;
end

% Convert cell to matrix in order to use the code of Kim, N.H. (2015)
Elem = cell2mat(Elem);

fprintf('Done!\n')
end