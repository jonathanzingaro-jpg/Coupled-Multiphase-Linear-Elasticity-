function Mesh = BuildMesh(Mesh)
% BUILDMESH Constructs the nodal locations and connectivity matrix for a
% structured mesh. The mesh is built for a model with spatial dimensions of
% nsd, length Lr, Lz, and ner, nez elements.  
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   nsd:    Number of spatial dimensions
%   L:      Length of domain (size nsd x 1)
%   ner:    Number of elements in each direction (size nsd x 1)
%   type:   Current assume Q4 elements 
% 
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   Mesh is a structure array with the following fields,
%       .nodeconn:      Array of nodal connectivity (size nn x 8)
%                       containing the indices of elements connected to 
%                       each node
%       .nsd:           Number of spatial dimensions
%       .ne:            Total number of elements in the mesh
%       .nne:           Vector of number of nodes per element (size nsd x 1)
%       .nn:            Total number of nodes 
%       .nodes:         Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)

%   Acknowledgments: Chris Ladubec, Matin Parchei Esfahani

% Print Message to Command Window
disp([num2str(toc),': Building Mesh...']);

% Mesh Properties from Mesh Input Structure
ner = Mesh.ner; nez = Mesh.nez; 
r_center = Mesh.r_center; % Radial Center of Wellbore
Lr = Mesh.Lr; Lz = Mesh.Lz;
nsd = Mesh.nsd; % No. of DOFs per node

if Mesh.type == 'Q4'

    % Mesh Properties
    ne = ner*nez;  % Total No. of Elements
    nnr = (ner+1); % number of nodes in each direction
    nnz = (nez+1); % number of nodes in each direction
    nn = nnr*nnz; % total number of nodes
    dLr = (Lr-r_center)/ner;% element size
    dLz = Lz/nez;% element size

    % Define Nodal Locations
    disp([num2str(toc),': Computing nodal locations...']);
    sr = r_center:dLr:Lr;
    sz = 0:dLz:Lz;
    nodes=zeros(nn,2);
    c=1;
    for i=1:nnz
        for j=1:nnr
            nodes(c,1)=sr(j);
            nodes(c,2)=sz(i);
            c=c+1;
        end
    end

    % Define Elemental Connectivity Matrix
    disp([num2str(toc),': Defining element connectivity...']);
    conn = zeros(ne,4);
    e=1;
    for i=1:nez
        for j=1:ner
            ii = (i-1)*nnr+(j-1)+1;
            conn(e,1)=ii;
            conn(e,2)=ii+1;
            conn(e,3)=ii+1+nnr;
            conn(e,4)=ii+nnr;
            e=e+1;
        end
    end

    disp([num2str(toc),': Saving mesh structure ...']);

    % Saving Mesh Properties
    Mesh.conn = conn; % Connectivity Matrix
    Mesh.nodes = nodes; % Spatial Coordinates of Nodes
    Mesh.nsd = nsd; % No. of DOFs per node
    Mesh.nne = 4; % No. of nodes per element
    Mesh.nn = nn; % Total no. of nodes
    Mesh.ne = ne; % No. of Total Elements in Mesh
    Mesh.eMat = ones(ne,1); % Assign Material Properties

elseif Mesh.type == 'Q9'

    % Mesh Properties
    ne = ner*nez;  % Total No. of Elements
    nnr = (2*ner+1); % number of nodes in each direction
    nnz = (2*nez+1); % number of nodes in each direction
    nn = nnr*nnz; % total number of nodes
    dLr = (Lr-r_center)/ner;% element size
    dLz = Lz/nez;% element size
    
    % Define Nodal Locations
    disp([num2str(toc),': Computing nodal locations...']);
    sr = r_center:dLr*0.5:Lr;
    sz = 0:dLz*0.5:Lz;
    nodes = zeros(nn,2);
    c = 1;
    for i=1:nnz
        for j=1:nnr
            nodes(c,1)= sr(j);
            nodes(c,2)= sz(i);
            c= c+1;
        end
    end

    % Define Elemental Connectivity Matrix
    disp([num2str(toc),': Defining element connectivity...']);
    conn = zeros(ne,9);
    e = 1;
    for i = 1:nez
        for j = 1:ner
            ii = 2*(i-1)*nnr + 2*(j-1)+1;
            conn(e,1) = ii;
            conn(e,2) = ii + 2;
            conn(e,3) = ii + 2 + 2*nnr(1);
            conn(e,4) = ii + 2*nnr(1);
            conn(e,5) = ii + 1;
            conn(e,6) = ii + 2 + nnr(1);
            conn(e,7) = ii + 1 + 2*nnr(1);
            conn(e,8) = ii + nnr(1);
            conn(e,9) = ii + 1 + nnr(1);
                
            e=e+1;
        end
    end

    % Saving Mesh Properties
    Mesh.conn = conn; % Connectivity Matrix
    Mesh.nodes = nodes; % Spatial Coordinates of Nodes
    Mesh.nsd = nsd; % No. of DOFs per node
    Mesh.nne = 9; % No. of nodes per element
    Mesh.nn = nn; % Total no. of nodes
    Mesh.ne = ne; % No. of Total Elements in Mesh
    Mesh.eMat = ones(ne,1); % Assign Material Properties

    disp([num2str(toc),': Saving mesh structure ...']);


end

disp([num2str(toc),': Done.']);

end
