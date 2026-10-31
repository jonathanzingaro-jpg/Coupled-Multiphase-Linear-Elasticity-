%% Disclosure 
%   Code runs and inputs parameters for a 1D Vertically Averaged Multi-
%   Phase Flow using SU Stabilization for a Single Injection Well
%   Created by Giovanni Zingaro for PhD Thesis

%% Clear Variables 
    % clear, clc, close all; % Clear the workspace, command window, and close all figures
    % format compact;       % Set the display format to compact
    % tic;                  % Start a timer

%% Change Current File Path

    % Current directory
    curDir = pwd;

    % Add Function Files to Current Directory
    FuncDir = fullfile(curDir, 'Functions');
    addpath(genpath(FuncDir));

    % Description
    description = '1D Vertically Averaged Multi-Phase Flow';

    % Output File Location
    data_dir = fullfile(curDir,'Data');

    % Output Messages 
    progress_on = 0;
    output = 0;

%% Material Inputs

    % Properties of CO2
    Material.muC = 1; % Viscosity of CO2 [Ns/m^2]
    Material.rhoC = 1; % Density of CO2 [kg/m^3]
    Material.krelC = 1; % Relative Permability of CO2 [-]

    % Properties of Brine
    Material.muB = 1; % Viscosity of Brine [Ns/m^2]
    Material.rhoB = 1; % Density of Brine [kg/m^3]
    Material.sB = 0; % Saturation of Brine [-]

    % Reservoir Properties 
    Material.kx = 1; % Permability [m^2]
    Material.phi = 1; % Porosity [-]

    % Other Properties 
    Material.g = 1; % Gravitational Constant [m/s^2]
    H = 1; % Height of Reservoir [m]
    Material.H = H;

%% Mesh Properties 
    Mesh.ner = N;       % No. of element in r direction
    Mesh.rw = 0;        % Radius of Wellbore [m] 
    Mesh.nsd = 1;         % No. of DOFs of Problem
    Mesh.L = 1;      % Length of Mesh
    Mesh.type = etype;     % Type of Elements
    Mesh.problemtype = 2; % 1 -> 2 Field Problem, 2 -> 1 Field Problem

    % Construct Mesh
    Mesh.FileType = 'MATLAB'; % File Types are (a) MATLAB, (b) Hypermesh -> MATLAB only supports Q4

    % Run Mesh Generator
    run('MeshGeneration.m')

%% Boundary Conditions

    % Reservoir Thickness 
    BC.H = H*ones(Mesh.nDOF,1);

    % Initial Pressure p0
    p0 = 0; 

    % Pressure Equation 

    % Dirichlet Boundary Condition (Essential)
    % ----------------------------------------

        % Column Vector of Prescribed Pressure DOF  
        BC.fixed_pres_dof = [Mesh.left_dof; Mesh.right_dof];

        % Prescribed Pressure for each DOF [u1; u2; ...] [Pa]
        BC.fix_pres_value = [];
       
    % Neumann Boundary Condition (Natural)
    % ---------------------------------------
        % Column Vector of Prescribed Traction DOF
        BC.flux_dof = [];

        % Magnitude of Prescribed Flux [m/s]
        BC.flux_dof_value = [];
 
    % Saturation Equation 
    
    % Dirichlet Boundary Condition (Essential)
    % ----------------------------------------

        % Column Vector of Prescribed Saturation DOF (hB)  
        BC.fixed_sat_dof = [Mesh.left_dof; Mesh.right_dof];

        % Prescribed Pressure for each DOF [u1; u2; ...] [Pa]
        BC.fixed_sat_value = [];  
        BC.fixed_sat_func  = @(r,t) 0.5*sin(pi.*r).*cos(pi.*t);

    % Neumann Boundary Condition (Natural)
    % ---------------------------------------
        % Column Vector of Prescribed Traction DOF
        BC.flux_dof = [];

        % Magnitude of Prescribed Flux [m/s]
        BC.flux_dof_value = []; 

    % Initial Condition
    t0 = 0;
    BC.IC_hb = (0.5.*sin(pi.*(Mesh.r)).*cos(pi.*t0)); % Initial Height of Brine Phase

    % Boundary Flux Functions
    BC.qc = @(r,t) 0;
    BC.qb = @(r,t) 0;
    
    % Source Function -> Pressure Equation
    BC.s = @(r,t) - ((r*pi^2*cos(pi*r)^2*cos(pi*t)^2)/2 - (r*pi^2*cos(pi*t)^2*sin(pi*r)^2)/2 + (pi*cos(pi*r)*cos(pi*t)^2*sin(pi*r))/2)/r - (pi*sin(pi*r)*sin(pi*t))/2;
   
    % Manufactured Solutions
    Control.Manufactured = 1;                                                    % ON: 1, OFF: 0
    Control.Manpbot = @(r,t) sin(pi.*r).*cos(pi.*t);
    Control.Manhb = @(r,t) (0.5.*sin(pi.*r).*cos(pi.*t));                     % Exact Solution of hb(r,t)
    Control.plot_toggle = 1;                                                     % ON: 1, OFF: 0

%% Computational Control
    % quadrature order
    Control.qo = 40;

    Q0 = 0;

    % Target Courant Number
    Control.Courant = 1;

    % Stabilization Method and Scale
    Control.Stab = 1; % 1 -> SU, 2-> GLS, 3 -> SUPG
    Control.Stab_Scale = 1; 

    % Parallel Inversion, use parallel processing to invert the matrix
    % Usually more efficient at 2e5 dofs
    Control.parallel = 1;

    % Control for Time Discritization
    Control.StartTime = 0;    % Start Time
    Control.EndTime   = 0.25; % End Time [s]
    Control.theta     = 1.0;  % Time Integration Parameter

    % Control for Saving Times
    Control.SaveTime = 10;
    Control.Save = 0;
    
    % Toggle for Coupled Mechanical Problem 
    Control.Mech = 0; % ON: 1, OFF: 0s
    mech_path = '\\galafrey\D$\U\gzingaro\Documents\Axisymmetric Elasticity Model';
      
%% Get Directory and Run Main File
run('main_Man.m')
