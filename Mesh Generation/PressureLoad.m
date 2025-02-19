% PressureLoad ------------------------------------------------------------
% Article: Fontes, V.O., Leitão, A.X., & Pereira, A. (2023). 
%          HyperSym: a Matlab code for symbolic differentiation of 
%          ready-to-use tensors in hyperelastic finite element analysis, 
%          Computer Applications in Engineering Education.
%          DOI: https://doi.org/xx.xxxx/xxxxxx-xxx-xxxx-x
% -------------------------------------------------------------------------
function EXTFORCE = PressureLoad(Nodes,LE,eRef,p,DISPTD)
% Compute pressure in the deformed configuration.
% The theoretical value for total nodal force is 2/4*pi*rt*0.775*150.
% To check this result use, sum(hypot(F(:,1),F(:,2))) or
% sum(vecnorm(F,2,2)).

if nargin<5, DISPTD = zeros(3*size(Nodes,1),1); end

Elem = LE(1:eRef(1):end,:);  el = 1;
F = zeros(size(Nodes));      LoadedNodes = zeros(4*size(Elem,1),1);

for IE = 1:size(LE,1)
    a = LE(IE,1); b = LE(IE,2); c = LE(IE,3); d = LE(IE,4);
    e = LE(IE,5); f = LE(IE,6); g = LE(IE,7); h = LE(IE,8);
    
    Faces = [...
        a b c d;
        e f g h;
        a b f e;
        d c g h;
        a d h e;   % <-- desired face
        b c g f];

    while Elem(el) == LE(IE)
        nd = Faces(5,:);          % nodes in the pressurized face
        nnd = length(nd);         % number of nodes in the face

        % Nodal coord. where the pressure is applied
        U = (reshape(DISPTD(:,end),3,[]))';
        Pt = Nodes(nd,:) + U(nd,:);
        
        % Vectors in the plane
        v1 = Pt(2,:) - Pt(1,:);
        v2 = Pt(3,:) - Pt(1,:);
        v3 = Pt(4,:) - Pt(1,:);
        
        if det([v1; v2; v3]) <= 1E-6
            da = cross(v1,v3);    % outward normal vector, deformed config.
            
            force = p*da/nnd;               % nodal force vector
            force = repmat(force,nnd,1);    % force distributed over nodes

            F(nd,:) = F(nd,:) + force;      % nodal force vector

            LoadedNodes(4*el-3:4*el,:) = nd;% pressurized nodes
        else
            error('Wrong reference point in the plane!')
        end
        
        el = el+1;
    end
end

% Eliminating unpressurized nodes (no pressure in z-direction)
F = F(any(F' ~= 0),1:2);

PressNodes = unique(LoadedNodes);

% Get the nodal force matrix to be used in FEA
EXTFORCE = [PressNodes   ones(length(PressNodes),1) F(:,1);...
            PressNodes 2*ones(length(PressNodes),1) F(:,2)];
end