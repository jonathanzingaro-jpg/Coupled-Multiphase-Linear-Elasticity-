function Mesh = NodeSets(Mesh)
%   NODESETS define sets of nodes and DOFs in the domain
%   Mesh = NODESETS(Mesh) updates the structure array containing mesh
%   information with relevant node sets
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:
%       .nsd        number of spatial dimensions
%       .nodes          array of nodal spatial locations for undeformed mesh 
%                   (size nn x nsd in which nn is the number of nodes)
% 
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   The function returns the Mesh structure array with new fields,
%       .left_nodes     Nodes on the left edge of the domain
%       .left_dof       DOFs on the left edge of the domain
%       .right_nodes    Nodes on the right edge of the domain
%       .right_dof      DOFs on the right edge of the domain
%       .xdofs          DOFs in the x-direction
%       .ydofs          DOFs in the y-direction
%       .zdofs          DOFs in the z-direction
%   Two-dimensional meshes also contain the new fields,
%       .top_nodes      Nodes on the top edge of the domain
%       .top_dof        DOFs on the top edge of the domain
%       .top_dofx       DOFs on the top boundary in the x-direction
%       .top_dofy       DOFs on the top boundary in the y-direction
%       .bottom_nodes   Nodes on the bottom edge of the domain
%       .bottom_dof     DOFs on the bottom edge of the domain
%       .bottom_dofx    DOFs on the bottom boundary in the x-direction
%       .bottom_dofy    DOFs on the bottom boundary in the y-direction
%       .left_dofx      DOFs on the left boundary in the x-direction
%       .left_dofy      DOFs on the left boundary in the y-direction
%       .right_dofx     DOFs on the right boundary in the x-direction
%       .right_dofy     DOFs on the right boundary in the y-direction
%   Three-dimensional meshes also contain the new fields, 
%       .near_nodes     nodes on the nearest face of the domain
%       .near_dof       DOFs on the nearest face of the domain
%       .near_dofx      DOFs on the near face in the x-direction
%       .near_dofy      DOFs on the near face in the y-direction
%       .near_dofz      DOFs on the near face in the z-direction
%       .far_nodes      Nodes on the farthest face of the domain
%       .far_dof        DOFs on the farthest face of the domain
%       .far_dofx       DOFs on the far face in the x-direction
%       .far_dofy       DOFs on the far face in the y-direction
%       .far_dofz       DOFs on the far face in the z-direction
%       .left_dofz      DOFs on the left face in the z-direction
%       .right_dofz     DOFs on the right face in the z-direction
%       .top_dofz       DOFs on the top face in the z-direction
%       .bottom_dofz    DOFs on the bottom face in the z-direction

% Acknowledgemnents: Chris Ladubec, Matin Parchei-Esfahani; Modified by
% Jonathan Zingaro


disp([num2str(toc),': Defining Node Sets ...']);

% Nodes and dofs on the left face of the domain
% (left and right are defined along the x-direction)
Mesh.left_nodes = find(Mesh.nodes(:,1)==min(Mesh.nodes(:,1)));
Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:);

% Nodes and dofs on the right face of the domain
Mesh.right_nodes = find(Mesh.nodes(:,1)==max(Mesh.nodes(:,1)));
Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);

% Nodes and dofs on the top face of the domain
% (top and bottom are defined along the y-direction)
Mesh.top_nodes = find(Mesh.nodes(:,2)==max(Mesh.nodes(:,2)));
Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);

% Nodes and dofs on the bottom face of the domain
Mesh.bottom_nodes = find(Mesh.nodes(:,2)==min(Mesh.nodes(:,2)));
Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);

Mesh.xdofs = 1:Mesh.nDOF;
Mesh.ydofs = [];
Mesh.zdofs = [];

Mesh.left_dofx = Mesh.left_dof(:,1);
Mesh.left_dofy = Mesh.left_dof(:,2);

Mesh.right_dofx = Mesh.right_dof(:,1);
Mesh.right_dofy = Mesh.right_dof(:,2);

Mesh.top_dofx = Mesh.top_dof(:,1);
Mesh.top_dofy = Mesh.top_dof(:,2);

Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
Mesh.bottom_dofy = Mesh.bottom_dof(:,2);

% DOFs in the x- and y- directions
Mesh.xdofs = 1:2:Mesh.nDOF;
Mesh.ydofs = 2:2:Mesh.nDOF;
Mesh.zdofs = [];

Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);
Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);
Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);
Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);

disp([num2str(toc),': Done']);

end













