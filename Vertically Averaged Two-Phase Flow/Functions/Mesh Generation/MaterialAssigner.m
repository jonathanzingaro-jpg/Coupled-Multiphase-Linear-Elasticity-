function [Material,Mesh] = MaterialAssigner(Material, Mesh)
% MaterialAssigner 
%   [Material, Mesh] = MeshMaterialAssigner(Material, Mesh) is used to assign the 
%   material parameters based on inputs in the RunSim file based on spatial
%   coordinates for layered media. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   nMat:       Material type
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .nsd:   Number of spatial dimensions
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .Dtype: constitutive law for two-dimensional elasticity
%               .E: Modulus of elasticity
%               .nu: Poisson's ratio
%   --------------------------------------------------------------------
%

% Created by: Jonathan Zingaro 

Matlim = Material.z; % Limits of Material Layers in Z-Direction
nMat = zeros(Mesh.ne,1);

% Define Material Types as a Loop 
for i = 1:Mesh.ne 
    coord = Mesh.nodes(Mesh.conn(i,:)',:); % Extract Coordinates of Nodes of Element i
    zmax = max(coord(:,2)); % Maximum Z-Coordinate in Element
    Mat_log = zeros(length(Matlim),1); % Define Vector to Store Logicals for Material Types

    % Loop over Number of Materials
    for j = 1:height(Matlim)
        if zmax < Matlim(j)
            Mat_log(j,1) = 1; % Assign Logical
        else
        end
    end
    nMat(i,1) = find(Mat_log==1, 1, 'last' ); % Assign Material Type

end

Material.nMat = nMat; % Store Material Types within Structure

end
