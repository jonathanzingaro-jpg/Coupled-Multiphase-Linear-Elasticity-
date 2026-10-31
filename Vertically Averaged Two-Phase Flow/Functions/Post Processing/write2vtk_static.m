function write2vtk_static(description, vtk_dir, Mesh, Control, BC, d, strain, stress, Fint, Fext, d_exact, e_exact, s_exact, Pressure)
%WRITE2VTK_STATIC Exports results to VTK file 
%   WRITE2VTK_STATIC(Mesh, Control, fixedDOF, d, strain, stress, Fint, Fext) 
%   produces vtk files that can be opened in Paraview to visualize
%   the simulation results. The mesh, boundary conditions, and 
%   simulation parameters are described in the structure arrays
%   Mesh, fixedDOF, and Control, respectively. The simulation results are 
%   the vectors d, Fint, and Fext of size ndof x 1, containing the 
%   displacement, internal forces, and external forces, respectively. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:       Structure array with the following fields, 
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .conn:  Array of element connectivity (size ne x nne)
%               .xdofs: DOFs in the x-direction
%               .ydofs: DOFs in the y-direction
%               .zdofs: DOFs in the z-direction
%               .nn:    Total number of nodes 
%               .nsd:   Number of spatial dimensions
%               .nDOF:  Total number of DOFs
%               .ne:    Total number of elements in the mesh
%   Control:    Structure array with the following fields,
%               .config_name:   Name of the configuration file
%               .vtk_dir:       Directory where VTK files are stored
%               .stress_calc:   Calculation of values for discontinous variables
%                               ('none', 'nodal', 'center')
%   fixedDOF:   Row vector containing fixed degrees of freedom
%   d:          Column vector of displacements (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   strain:     if Control.stress_calc = 'nodal': Matrix of nodal strains 
%                           (size dim x nn in which dim is 1 for 1D 
%                           elements, 3 for 2D elements, and 6 for 3D 
%                           elements and nn is the number of nodes)
%               if Control.stress_calc = 'center': Matrix of element 
%                           strains (calculated at the center of each 
%                           element); size dim x ne in which ne is the 
%                           number of elements
%   stress:     Matrix of stresses (same size as strain matrix)
%   Fint:       Column vector of internal forces (size ndof x 1)
%   Fext:       Column vector of external forces (size ndof x 1)

%% Define variables
    filename = fullfile(vtk_dir, [description '.vtk']); 

    R = Fext - Fint;

%% Declare 
fixedDOF = BC.fixed;

%% Data
    nodedata(1).name = 'nID';
    nodedata(1).data = (1:Mesh.nn)';
    nodedata(1).type = 'int';

    elementdata(1).name = 'eID';
    elementdata(1).data = (1:Mesh.ne)';
    elementdata(1).type = 'int';

    nodedata(end+1).name = 'dof';
    nodedata(end).data = Mesh.DOF;
    nodedata(end).type = 'int';

    fixed = zeros(Mesh.nDOF,1);
    fixed(fixedDOF) = 1;
    nodedata(end+1).name = 'fixed';
    nodedata(end).data = fixed(Mesh.DOF);
    nodedata(end).type = 'int';

    % Store Displacement, Stress, Strain, Fext and Fint 
    nodedata(end+1).name = 'U';
    nodedata(end).data = [d(Mesh.DOF), zeros(size(Mesh.xdofs'))];
    nodedata(end).type = 'float';

    nodedata(end+1).name = 'Fext';
    nodedata(end).data = [Fext(Mesh.DOF), zeros(size(Mesh.xdofs'))];
    nodedata(end).type = 'float';

    nodedata(end+1).name = 'Fint';
    nodedata(end).data = [Fint(Mesh.DOF), zeros(size(Mesh.xdofs'))];
    nodedata(end).type = 'float';

    nodedata(end+1).name = 'R';
    nodedata(end).data = [R(Mesh.DOF), zeros(size(Mesh.xdofs'))];
    nodedata(end).type = 'float';


    if strcmp(Control.StressCalc,'Nodal')
        nodedata(end+1).name = 'err';
        nodedata(end).data = strain(1,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'ezz';
        nodedata(end).data = strain(2,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'ett';
        nodedata(end).data = strain(3,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'erz';
        nodedata(end).data = strain(4,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Srr';
        nodedata(end).data = stress(1,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Szz';
        nodedata(end).data = stress(2,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Stt';
        nodedata(end).data = stress(3,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Srz';
        nodedata(end).data = stress(4,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'p';
        nodedata(end).data = Pressure.p(Mesh.nodes(:,1),Mesh.nodes(:,2));
        nodedata(end).type = 'float';

        if Control.Manufactured == 1
            d_L2 = ((d(Mesh.DOF)-d_exact(Mesh.DOF)).^2);
            nodedata(end+1).name = 'L2-Ur';
            nodedata(end).data = d_L2(:,1);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'L2-Uz';
            nodedata(end).data = d_L2(:,2);
            nodedata(end).type = 'float';

            err_L2 = (strain(1,:)- e_exact(1,:)).^2;
            nodedata(end+1).name = 'L2-err';
            nodedata(end).data = err_L2';
            nodedata(end).type = 'float';

            ezz_L2 = (strain(2,:)- e_exact(2,:)).^2;
            nodedata(end+1).name = 'L2-ezz';
            nodedata(end).data = ezz_L2';
            nodedata(end).type = 'float';
            
            ett_L2 = (strain(3,:)- e_exact(3,:)).^2;
            nodedata(end+1).name = 'L2-ett';
            nodedata(end).data = ett_L2';
            nodedata(end).type = 'float';
      
            erz_L2 = (strain(4,:)- e_exact(4,:)).^2;
            nodedata(end+1).name = 'L2-erz';
            nodedata(end).data = erz_L2';
            nodedata(end).type = 'float';

            srr_L2 = (stress(1,:)- s_exact(1,:)).^2;
            nodedata(end+1).name = 'L2-srr';
            nodedata(end).data = srr_L2';
            nodedata(end).type = 'float';

            szz_L2 = (stress(2,:)- s_exact(2,:)).^2;
            nodedata(end+1).name = 'L2-szz';
            nodedata(end).data = szz_L2';
            nodedata(end).type = 'float';
            
            stt_L2 = (stress(3,:)- s_exact(3,:)).^2;
            nodedata(end+1).name = 'L2-stt';
            nodedata(end).data = stt_L2';
            nodedata(end).type = 'float';
      
            srz_L2 = (stress(4,:)- s_exact(4,:)).^2;
            nodedata(end+1).name = 'L2-srz';
            nodedata(end).data = srz_L2';
            nodedata(end).type = 'float';

        else

        end

    elseif strcmp(Control.stress_calc, 'center')

        nodedata(end+1).name = 'err';
        nodedata(end).data = strain(1,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'ezz';
        nodedata(end).data = strain(2,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'ett';
        nodedata(end).data = strain(3,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'erz';
        nodedata(end).data = strain(4,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Srr';
        nodedata(end).data = stress(1,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Szz';
        nodedata(end).data = stress(2,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Stt';
        nodedata(end).data = stress(3,:)';
        nodedata(end).type = 'float';

        nodedata(end+1).name = 'Srz';
        nodedata(end).data = stress(4,:)';
        nodedata(end).type = 'float';
    end

%% Write to vtk
    WriteMesh2VTK(filename, description, Mesh.nodes, ...
                    Mesh.conn, nodedata, elementdata);

end
