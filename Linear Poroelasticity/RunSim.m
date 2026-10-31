%% Disclosure 
%   Code runs and inputs parameters for an axisymmetric elasticity code
%   where pore pressure is defined from the 1D vertically averaged
%   multiphase fluid code  
%   Created by Giovanni Zingaro for PhD Thesis

%% Clear Variables 
    % clear, clc, close all
    % format compact 
    % tic; 

%% Change Current File Path

    % % Fit Pressure Build-up
    fileID = fopen(fullfile(filename,'pres.txt'),'r');      
    s = textscan(fileID, '%s', 'delimiter', '\n');
    t = s{1}{2,1};
    p = s{1}(6:end); % Extract Pressure Build Values
    p = str2double(p); % Convert to Double
    fclose(fileID); % Close File for Stability 


    % % Create VTK Name
    pattern = 'Time [days]: ';
    index = strfind(t, pattern);
    startIdx = index + length(pattern);
    extractedSubstring = num2str(floor(str2double(t(startIdx:end))));

    VTKname = ['WellboreVTK' extractedSubstring];

    [p_function] = fitPressure1D(p, disp_fig, r, 'cubicinterp'); % Type:'cubicinterp',

    % Current directory
    curDir = pwd;

    % Add Function Files to Current Directory
    FuncDir = fullfile(curDir, 'Functions');
    addpath(genpath(FuncDir));
 
    % VTK Directory
    VTKfilename = VTKfilename_path;
    description = [VTKname '.vtk'];

    % Data File Path 
    data_dir = filename;

%% Material Inputs

    % % Material Properties - Pre-Cambrian (Resevoir) (ID: 1 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(1).E = 100e9;     % Young's Modulus [Pa]
    % Material.Prop(1).nu = 0.25;     % Poisson's Ratio [-]
    % Material.Prop(1).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Cambrian (Resevoir) (ID: 2 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(2).E = 12e9;      % Young's Modulus [Pa]
    % Material.Prop(2).nu = 0.38;     % Poisson's Ratio [-]
    % Material.Prop(2).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Trenton and Black River Groups (ID: 3 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(3).E = 13.8e9;    % Young's Modulus [Pa]
    % Material.Prop(3).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(3).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Georgian Bay Formation (ID: 4 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(4).E = 9e9;       % Young's Modulus [Pa]
    % Material.Prop(4).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(4).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Queenston Formation (ID: 5 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(5).E = 15e9;      % Young's Modulus [Pa]
    % Material.Prop(5).nu = 0.4;      % Poisson's Ratio [-]
    % Material.Prop(5).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Clinton and Medina Groups (ID: 6 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(6).E = 13.8e9;    % Young's Modulus [Pa]
    % Material.Prop(6).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(6).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Lockport Group (ID: 7 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(7).E = 35e9;      % Young's Modulus [Pa]
    % Material.Prop(7).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(7).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Salina Group (ID: 8 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(8).E = 20e9;      % Young's Modulus [Pa]
    % Material.Prop(8).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(8).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Bass Island Formation (ID: 9 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(9).E = 4e9;       % Young's Modulus [Pa]
    % Material.Prop(9).nu = 0.3;      % Poisson's Ratio [-]
    % Material.Prop(9).b = [0; 0];    % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Devonian Age Groups (ID: 10 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(10).E = 7.7e9;    % Young's Modulus [Pa]
    % Material.Prop(10).nu = 0.18;    % Poisson's Ratio [-]
    % Material.Prop(10).b = [0; 0];   % Body Force Density [br; bz] [N/m^3]
    % 
    % % Material Properties - Pleistoncene Age Groups (ID: 11 from Hypermesh)
    % % -----------------------------------------------------
    % Material.Prop(11).E = 10e9;     % Young's Modulus [Pa]
    % Material.Prop(11).nu = 0.2;     % Poisson's Ratio [-]
    % Material.Prop(11).b = [0; 0];   % Body Force Density [br; bz] [N/m^3]
    % 
    % % Reservoir Properties 
    % Mesh.Hmax = 135;                % Maximum Height of Reservoir
    % Mesh.Hmin = 100;                % Minimum Height of Reservoir
    % 
    % % Constitutive Law -> Axisymmetric, Plane Strain, Plane Stress
    % Material.Dtype = 'Axisymmetric'; 

% Mesh Properties 

    % Mesh.nsd = 2; % No. of DOFs of Problem
    % Mesh.type = 'Q4'; % Supports T3, Q4 and Q8
    % 
    % % Construct Mesh
    problemtype = 1; 
    % Mesh.FileType = 'Hypermesh'; % File Types are (a) MATLAB, (b) Hypermesh -> MATLAB only supports Q4
    % Mesh.FileName = fullfile(curDir,'/Mesh/Subsurface - 100 m - 25km - Q4 - Variable Measurement Points.fem');
    % 
    % BC.Hypermesh_Measurement = 1;
    % 
    % % Assembly Mesh
    % run('MeshGeneration.m')

% Boundary Conditions
    
    % % Dirichlet Boundary Condition (Essential)
    % % ----------------------------------------
    % 
    %     % column vector of prescribed displacement dof  
    %     BC.fix_disp_dof = BC.fixed;
    % 
    %     % prescribed displacement for each dof [u1; u2; ...] [m]
    %     BC.fix_disp_value = zeros(height(BC.fixed),1);  
    % 
    % % Neumann Boundary Condition (Natural)
    % % ---------------------------------------
    %     % column vector of prescribed traction dofs
    %     BC.traction_force_dof = [];
    % 
    %     % magnitude of prescribed tractions [N]
    %     BC.traction_force_dof_value = [];
    % 
    %     % NOTE: this is slower than prescribing tractions at dofs
    %     % column vector of prescribed traction nodes 
    %     BC.traction_force_node = [];  
    % 
    %     % prescribed traction [t1x t1y;t2x t2y;...] [N]
    %     Fnode = 0;
    %     BC.traction_force_value = [];
    
    % Pore Pressure Term
    % --------------------------------------
        Pressure.p = p_function; % Pore Pressure Function from 1D Vertically Averaged Code


%% Computational Control
    % quadrature order
    Control.qo = 3;

    % parallel inversion
    % Use parallel processing to invert the matrix.
    % Usually more efficient at 2e5 dofs
    Control.parallel = 1;

    % Calculation of values for discontinuous variables 
    % (i.e. stress/strain)
    % 'none': calculated at each node for each element separately; 
    %           no output in vtk
    % 'nodal': averaged at each node for all elements attached to 
    %           the node; output as nodal values in vtk
    % 'center': calculated at the center of each element; output as 
    %           single value for each element in vtk
    % 'L2projection': Least squares projection of stress and strain,
    %           output as nodal values
    Control.StressCalc = 'Nodal';

    % Toggle for Manufactured Solution
    Control.Manufactured = 0; % ON = 1; OFF = 0;

    % Toggle for Generatoring Data
    Control.DataGenerator = 1; % ON = 1; OFF = 0;

    % Toggle for Hydrostatic Pore Pressure
    Control.Hydrostatic = 1;

    % Plot to VTK
    Control.vtk = 1;

    % Toggle for Stress and Strain
    Control.Stress_Toggle = 1;

    % Batch File Generation
    Control.Batch = 0;

%% Get Directory and Run Main File
run('main.m')
