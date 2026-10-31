%% Axisymmetric Elastic Material -> Main File
% Acknowledgements: Giovanni Zingaro

% % Identify Free and Fixed DOFs
% if Control.Batch == 0
%     BC.fixed = BC.fix_disp_dof;
%     BC.free = setdiff(Mesh.DOF, BC.fixed)';
% 
% %% Quadrature calculation
%     % (Done at beginning of code - this assumes all elements are the 
%     % same type and same order of quadrature for all integrals)
%     Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);
% 
% %% Compute System Matrices and Vectors
%     % Linear Elastic Siffness Matrix in Axisymmetric Coordinate System
%     disp([num2str(toc),': Assembling Linear Elastic Stiffness Matrix...']);
%     K = getK(Mesh,Quad,Material);
%     disp([num2str(toc),': Done']);
% end

    % Compute External Force Vector 
    disp([num2str(toc),': Assembling External Force Vector...']);
    Fext = getFext(Mesh, BC, Quad, Control, Material, Pressure);
    disp([num2str(toc),': Done']);    

%% Solve System of Equations 

    % Initialize Displacement Values 
    d = zeros(Mesh.nDOF,1);

    % Prescribe Essential Boundary Conditions
    d(BC.fix_disp_dof) = BC.fix_disp_value;

    % Linear Solver
    disp([num2str(toc),': Solving System of Equations...']);
    [d,fE] = LinearSolver1(K,Fext,d(BC.fix_disp_dof),BC.free,BC.fixed,Control.parallel);
    disp([num2str(toc),': Done']); 

    % Interal Force (Fint)
    Fint = K*d;

    % Compute Stress and Strain
    if Control.Stress_Toggle == 1
        disp([num2str(toc),': Computing Stress and Strain...']);
        [stress,strain] = getStrain(d, Mesh, Material, Control.StressCalc, Quad);
        disp([num2str(toc),': Done']);
    end

%% Compare Manufactured Solution

%% Compare Manufactured Solution
    % if Control.Manufactured == 1
    %     disp([num2str(toc),': Checking Manufactured Solution...']);
    %     [m_L2,m_e, d_exact, e_exact, s_exact] = ManufacturedSolutionCheck(d,stress, strain, Mesh, Control, Material);
    %     store(end+1,:) = [m_L2,m_e];
    %     disp([num2str(toc),': Done']); 
    % else
    % end

%% Write to VTK
    
    if Control.vtk == 1

        if Control.Manufactured == 0
            d_exact = 0; e_exact = 0; s_exact = 0;
        end

        disp([num2str(toc),': Writing to VTK...']);
        write2vtk_static(description, VTKfilename, Mesh, Control, ...
                        BC, d, strain, stress, Fint, Fext, d_exact, e_exact, s_exact, Pressure)
        disp([num2str(toc),': Done']); 
    end

%% Export Data to Textfile
    %if Control.DataGenerator == 1
        disp([num2str(toc),': Writing to Data to Textfile...']);
        writeData2txt_Disp(d, Mesh, data_dir, BC)
        disp([num2str(toc),': Done']);
    %end
