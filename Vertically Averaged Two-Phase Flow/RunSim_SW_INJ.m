%% Disclosure 
%   Code runs and inputs parameters for a 1D Vertically Averaged Multi-
%   Phase Flow using SU Stabilization for a Single Injection Well
%   Created by Giovanni Zingaro for PhD Thesis

%% Clear Variables 
    % clear, clc, close all; % Clear the workspace, command window, and close all figures
    % format compact;       % Set the display format to compact
    % tic;                  % Start a timer

%% Change Current File Path

    % Display parameters
    disp('Running simulation with parameters:');
    disp(['tup = ', num2str(tup)]);
    disp(['tend = ', num2str(tend)]);
    disp(['Q0 = ', num2str(Q0)]);
    disp(['dQup = ', num2str(dQup)]);
    disp(['k = ', num2str(k)]);
    disp(['phi = ', num2str(phi)]);
    disp(['iter = ', num2str(iter)]);

    % % Modify Values to Align with Units in Code 
    % dtQ_i = 91.25; % Number of Days before dQ
    % tup = tup * 60 * 60* 24 * 365;
    % tend = tend * 60 * 60 * 24 * 365 ;
    % dtQ = dtQ_i * 60 * 60 * 24;
    % Q0 = 0.0022 * (Q0).*(1e9)./(2 * 60 * 60 * 24 * 365);
    % dQup = 0.0022 * (dQup).*(1e9)./(2 * 60 * 60 * 24 * 365) ./ (tup);
    
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
    Material.muC = 0.062e-3;            % Viscosity of CO2 [Ns/m^2]
    Material.rhoC = 750;                % Density of CO2 [kg/m^3]

    % Properties of Brine
    Material.muB = 0.84e-3;             % Viscosity of Brine [Ns/m^2]
    Material.rhoB = 1045;               % Density of Brine [kg/m^3]
    Material.sB = 0.3;                  % Saturation of Brine [-]

    % Reservoir Properties 
    Material.kx = k * 9.869233e-16;     % Permability [m^2]
    Material.phi = phi;                 % Porosity [-]
    Material.krelC = 0.8;               % Relative Permeability of CO2 [-]

    % Other Properties 
    Material.g = 9.81;                  % Gravitational Constant [m/s^2]
    Mesh.H = 35;                        % Height of Reservoir [m]    Mesh.H = 35;                        % Height of Reservoir [m]
    Material.Bg = 0.0022;               % Change in Volume of CO2 from Surface to Reservoir 

%% Mesh Properties 
    Mesh.ner = 1000;                    % No. of element in r direction
    Mesh.rw = 0.2;                      % Radius of Wellbore [m] 
    Mesh.nsd = 1;                       % No. of DOFs of Problem
    Mesh.L = 25*1e3;                    % Length of Mesh
    Mesh.type = 'L2';                   % Type of Elements
    Mesh.problemtype = 2;               % 1 -> 2 Field Problem, 2 -> 1 Field Problem

    % Construct Mesh
    Mesh.FileType = 'MATLAB';           % File Types are (a) MATLAB

    % Run Mesh Generator
    run('MeshGeneration.m')

%% Boundary Conditions

    % Elevation Properties
    BC.H = Mesh.H*ones(Mesh.nDOF,1);

    % Pressure Equation 

    % Dirichlet Boundary Condition (Essential)
    % ----------------------------------------

        % Column Vector of Prescribed Pressure DOF  
        BC.fixed_pres_dof = [Mesh.right_dof];

        % Prescribed Pressure for each DOF [u1; u2; ...] [Pa]
        BC.fix_pres_value = 0; 
       
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
        BC.fixed_sat_dof = [Mesh.right_dof];

        % Prescribed Pressure for each DOF [u1; u2; ...] [Pa]
        BC.fixed_sat_value = [Mesh.H];  

    % Neumann Boundary Condition (Natural)
    % ---------------------------------------
        % Column Vector of Prescribed Traction DOF
        BC.flux_dof = [Mesh.left_dof];

        % Magnitude of Prescribed Flux [m/s]
        BC.flux_dof_value = []; 

    % Initial Condition
    BC.IC_hb = (Mesh.H)*ones(height(Mesh.r),1); % Initial Height of Brine Phase

    % Injection Nodes
    BC.qb = @(r,t) 0;

    % Injection Source Function
    BC.s = @(r,t) 0;

%% Computational Control
    % quadrature order
    Control.qo = 4;

    % Target Courant Number
    Control.Courant = 5;

    % Stabilization Method and Scale
    Control.Stab = 1; % 1 -> SU, 2-> GLS, 3 -> SUPG
    Control.Stab_Scale = 75; 

    % Parallel Inversion, use parallel processing to invert the matrix
    % Usually more efficient at 2e5 dofs
    Control.parallel = 1;

    % Toggle for Manufactured Solution
    Control.Manufactured = 0; %

    % Control for Time Discritization
    Control.StartTime = 0;    % Start Time
    Control.EndTime   = tend; % End Time [s]
    Control.theta     = 1.0;  % Time Integration Parameter

    % Control for Saving Times
    Control.Save = 1;
    Control.SaveTime = sort([[0,dtQ:dtQ:tend],tend]);

    % Toggle for Coupled Mechanical Problem 
    Control.Mech = 0; % ON: 1, OFF: 0s
    mech_path = '/Users/jonathanzingaro/Documents/PhD/Research/Axisymmetric Elasticity Model/Axisymmetric Elasticity Model';
      
%% Get Directory and Run Main File
run('main.m')
