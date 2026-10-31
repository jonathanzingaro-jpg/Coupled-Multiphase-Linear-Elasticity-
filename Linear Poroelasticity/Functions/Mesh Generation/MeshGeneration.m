% Axisymmetric Elastic Material -> Mesh Generation File
% Acknowledgements: Giovanni Zingaro

%% Mesh Assembly -> Identify Connectivity Matrix, Node Location

switch Mesh.FileType
    case 'MATLAB'
        Mesh = BuildMesh(Mesh);
        %Mesh = BuildMesh_structured(Mesh.nsd, ner, Mesh.type, problemtype, Mesh, 1, [Mesh.r_center, 0]);
        
    case 'Hypermesh'
        [Mesh, BC] = BuildMesh_Hypermesh(Mesh);
end

% Assign Global Nodal No.
Mesh = NodeDOFs(Mesh, problemtype);

% Group Node Sets for Boundary Conditions
Mesh = NodeSets(Mesh);


