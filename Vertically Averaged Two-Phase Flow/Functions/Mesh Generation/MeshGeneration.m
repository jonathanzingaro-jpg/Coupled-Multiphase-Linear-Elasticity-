% Axisymmetric Elastic Material -> Mesh Generation File
% Acknowledgements: Giovanni Zingaro

%% Mesh Assembly -> Identify Connectivity Matrix, Node Location

switch Mesh.FileType
    case 'MATLAB'
        Mesh = BuildMesh_structured(Mesh, progress_on);
        
    case 'Hypermesh'
        [Mesh, BC] = BuildMesh_Hypermesh(Mesh, progress_on);
end

% Assign Global Nodal No.
Mesh = NodeDOFs(Mesh);

% Group Node Sets for Boundary Conditions
Mesh = NodeSets(Mesh);


