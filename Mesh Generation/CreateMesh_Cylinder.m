% CreateMesh_Cylinder -----------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function [Nodes,Elem] = CreateMesh_Cylinder(Dim,BdBox,NElem)
% Generate mesh. It starts from a 2D model which is extruded.
% Then the element connectivities are performed.

NElemX = NElem(1); % number of elements along X-axis
NElemY = NElem(2); % number of elements along Y-axis
NElemZ = NElem(3); % number of elements along Z-axis
          
DeltaX = (BdBox(2) - BdBox(1))/NElemX;
DeltaZ = (BdBox(6) - BdBox(5))/NElemZ;

fprintf('Generating mesh...\n\n')

% 2D Mesh -----------------------------------------------------------------
% Node Assembling
NodesInPlane = zeros((NElemX + 1)*(NElemY + 1),2);
theta = 0;

for ndy = 1:(NElemY + 1)
    for ndx = 1:(NElemX + 1)
        NodesInPlane((NElemX + 1)*(ndy - 1) + ndx,:) = [...
            (BdBox(1) + (ndx - 1)*DeltaX)*cosd(theta) ...
            (BdBox(3) + (ndx - 1)*DeltaX)*sind(theta)];
    end
    
    theta = theta + Dim.alpha/NElemY;
end

% Element Assembling
ElemInPlane = cell(NElemX*NElemY,1);

for elx = 1:NElemX
    for ely = 1:NElemY
        ElemInPlane{NElemX*(ely - 1) + elx} = [...
            (NElemX + 1)*(ely - 1) + elx ...
            (NElemX + 1)*(ely - 1) + elx + 1 ...
            (NElemX + 1)*ely + elx + 1 ...
            (NElemX + 1)*ely + elx];
    end
end

% 3D Mesh -----------------------------------------------------------------
% Node Assembling
Nodes = zeros((NElemX + 1)*(NElemY + 1)*(NElemZ + 1),3);
Nodes(:,1:2) = repmat(NodesInPlane,NElemZ+1,1); % extrusion in z-direction

row = 0;     t = 0;     Ne = length(NodesInPlane);

while row <= NElemZ     % fill the extruded coordinates
    Nodes(row*Ne + 1:(row+1)*Ne,3) = t;
    t = t + DeltaZ;     row = row + 1;
end

% Element Assembling
NNodesInPlane = length(NodesInPlane);
Elem = cell(NElemX*NElemY*NElemZ,1);
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