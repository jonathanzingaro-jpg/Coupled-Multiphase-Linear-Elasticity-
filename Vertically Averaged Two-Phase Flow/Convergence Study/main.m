%% Axisymmetric Vertically Averaged Multi-Phase Flow Code -> Main File
% Acknowledgements: Giovanni Zingaro

%% Declare Function Pointers for Stabilization Method 
    [storagestabfile_name, diffusionstabfile_name] = setStabilizationMethod(Control);

%% Natural and Essential Boundary Conditions
    
    % Define Boundary Conditions for Saturation Equation
    BC.fixed_sat = BC.fixed_sat_dof;
    BC.free_sat = setdiff(Mesh.DOF, BC.fixed_sat)';

    % Define Boundary Conditions for Pressure Equation
    BC.fixed_pres = BC.fixed_pres_dof;
    BC.free_pres = setdiff(Mesh.DOF, BC.fixed_pres)';
    
%% Quadrature calculation
    % (Done at beginning of code - this assumes all elements are the 
    % same type and same order of quadrature for all integrals)
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);

%% Set Initial Condition for hb 
    h.hbn1 = BC.IC_hb;
    h.hcn1 = H*ones(Mesh.nn,1)-h.hbn1;
    h.hc = h.hcn1;
    h.hb = h.hbn1;

%% Initialize Time and Spatial Parameters
    t = Control.StartTime; % Current Time
    tend = Control.EndTime; % End Time
    theta = Control.theta; % Numerical Integration Parameter
    dr = Mesh.L/Mesh.ner; % Spatial Discretization
    Qt = Q0; % Initialize Injection Rate

    % Update Boundary Condition for Flux of CO2
    if Control.Manufactured == 0
        BC.qc = @(r,t) Q0;
    end

%% Storage Arrays
    tstore = [];
    Qstore = [];
    pwstore = [];

%% Assemble and Solve Pressure Equations for Initial Condition
    
    % Apply Essential BC for Manufactured Solution
    if Control.Manufactured == 1
        BC.fix_pres_value = Control.Manpbot((Mesh.r(BC.fixed_pres_dof)),t);
    end

    % Linear Diffusion/Stiffness Matrix for Pressure Equation
    Kp = getKp(Mesh, Quad, Material, h); % Diffusion/Stiffness Matrix

    % Boundary Flux Vector for Pressure Equation
    Fp1 = getFp1(Mesh, t, BC); 

    % Injection Flux Vector for Pressure Equation
    Fp2 = getFp2(Mesh, Quad, BC, t);

    % Buonyancy Vector for Pressure Equation
    Fp3 = getFp3(Mesh, Quad, h.hc, Material);

    % Elevation Flux Vector 
    Fp4 = getFp4(Mesh, Quad, h, Material);

    % Sum the External Flux Vector
    Fp = Fp1 + Fp2 - Fp3 + Fp4;%zeros(Mesh.nn,1);

    % Solve the System of Equations for the Pressure Equation
    [pbot,~] = LinearSolver1(Kp,Fp,BC.fix_pres_value,BC.free_pres,BC.fixed_pres,Control.parallel);

    % Save Wellbore Pressure
    pwstore(end+1,1) = pbot(1); 

%% Declare Timestep for First Step
    Control.TimeStep = 1e-10;
    dt = Control.TimeStep;
    counter = 1; % Timestep Counter 

    %% Save Intial Conditions

    % Save Counter 
    save_counter = 1;

    % Save Initial Conditions
    if Control.Save == 1
        writeData2txt(pbot, h.hb, Mesh, Material, data_dir, iter, Qt, t, counter, BC);
    end

    % Run Mechancial Code
    if Control.Mech == 1
        RunSimMech(Mesh.r,pbot,t,data_dir,iter,mech_path)

        % Update Function Paths
        curDir = pwd;

        % Add Function Files to Current Directory
        FuncDir = fullfile(curDir, 'Functions');
        addpath(genpath(FuncDir));

    end

    if (Control.Save) == 1

        % Next Saving Time
        save_counter = save_counter + 1;
        tsave = Control.SaveTime(save_counter);

    end

    % Initialize Injection Vector for Saturation Equation
    Fgn1 = zeros(Mesh.nDOF, 1);

%% Precompute Storage Matrix for Saturation Equation 
    Cstab = getCsu(Mesh, Quad, pbot, Material, Control);

%% Timestep Loop
    while t<=tend

        %% Timestep
        % If tend is execeed reduce timestep size

        %% Assemble and Solve Pressure Equations

        % % Apply Essential BC for Manufactured Solution
        % if Control.Manufactured == 1
        %     BC.fix_pres_value = Control.Manpbot((Mesh.r(BC.fixed_pres_dof)),t);
        % end
        % 
        % % Update Boundary Condition for Flux of CO2
        % if Control.Manufactured == 0
        %     BC.qc = @(r,t) Qt;
        % end
        % 
        % % Linear Diffusion/Stiffness Matrix for Pressure Equation
        % Kp = getKp(Mesh, Quad, Material, h); % Diffusion/Stiffness Matrix
        % 
        % % Boundary Flux Vector for Pressure Equation
        % Fp1 = getFp1(Mesh, t, BC); 
        % 
        % % Injection Flux Vector for Pressure Equation
        % Fp2 = getFp2(Mesh, Quad, BC, t);
        % 
        % % Buonyancy Vector for Pressure Equation
        % Fp3 = getFp3(Mesh, Quad, h.hc, Material);
        % 
        % % Elevation Flux Vector
        % Fp4 = getFp4(Mesh, Quad, h, Material);
        % 
        % % Sum the External Flux Vector
        % Fp = Fp1 + Fp2 - Fp3 + Fp4;
        % 
        % % Solve the System of Equations for the Pressure Equation
        % [pbot,~] = LinearSolver1(Kp,Fp,BC.fix_pres_value,BC.free_pres,BC.fixed_pres,Control.parallel);

        % Exact Pressure Build-up
        pbot = sin(pi.*(Mesh.r)).*cos(pi.*t);

        % Compute Minimum Timestep
        dt = ComputeTimeStep(Mesh, Quad, Material, pbot, Control, BC);
        Control.TimeStep = dt;

        if t == tend
            break
        elseif t+dt > tend
            dt = tend-t;
            % elseif (t+dt >= tsave)
            %     dt = tsave-t;
            %
            %     % Modify Flow Rate
            %     if t <= tup
            %         Qt = Qt + dQup * (Control.SaveTime(save_counter+1) - Control.SaveTime(save_counter));
            %
            %     elseif (t > (tup+tconst))
            %         Qt = max(Qt - dQdown * (Control.SaveTime(save_counter+1) - Control.SaveTime(save_counter)),0);
            %     end
            %
            %     % % Store Time and Injection Rate
            %     % tstore(end+1,1) = t;
            %     % Qstore(end+1,1) = Qt;

        else
        end

        t = t + dt; % Update Current Time

        % Exact Pressure Build-up
        pbot = sin(pi.*(Mesh.r)).*cos(pi.*t);

        % Assemble Fully Stabilized Advection Reaction Matrix (K)
        [Ks, dtn1] = getKsu(Mesh, Quad, Material, pbot, Control);

        % Update Timestep
        % Control.TimeStep = dtn1;
        % dt = dtn1;

        % Assemble Stabilized Storage Matrix (Cstab)
        C = Cstab; % Fully Stabilized Storage Matrix

        % Declare Flux Vector for Saturation Equation
        Fg = getFg1(Mesh, Quad, BC, t);

        %% Solve Saturation Equation

        % Compute Essential BC
        if Control.Manufactured == 1
            BC.fixed_sat_value = Control.Manhb((Mesh.r(BC.fixed_sat_dof)),t);
        end

        % Solve the Transient Problem
        Keff = C+dt*theta*Ks;
        Feff = (C-dt*(1-theta)*Ks)*(h.hbn1)+dt*(theta*Fg+(1-theta)*Fgn1);

        % Update Injection Vector 
        Fgn1 = Fg;

        % Solve the System of Equations for the Saturation Equation
        [hb,fE] = LinearSolver1(Keff,Feff,BC.fixed_sat_value,BC.free_sat,BC.fixed_sat,Control.parallel);

        % % Slope Limiter 
        % negativeIndex = find(hb < 0, 1, 'first');

        % % % Check if a negative value was found
        % if ~isempty(negativeIndex)
        % 
        %     % Set all values from the first negative index to the end of the array to 0
        %     hb(1:negativeIndex) = 0;
        % 
        % end

        % Update hb and hc
        h.hbn1 =  hb;
        h.hcn1 = H*ones(Mesh.nn,1) - hb;

        % % Save to Textfile
        % if Control.Save == 1
        %     if t == tsave
        % 
        %         % Save data to text file
        %         writeData2txt(pbot, hb, Mesh, Material, data_dir, iter, Qt, t, counter, BC);
        % 
        %         % Run Mechancial Code
        %         if Control.Mech == 1
        %             RunSimMech(Mesh.r,pbot,t,data_dir,iter,mech_path)
        % 
        %             % Update Function Paths
        %             curDir = pwd;
        % 
        %             % Add Function Files to Current Directory
        %             FuncDir = fullfile(curDir, 'Functions');
        %             addpath(genpath(FuncDir));
        % 
        %         end
        % 
        % 
        %         if t == Control.SaveTime(end)
        %         else
        %             save_counter = save_counter + 1;
        %             tsave = Control.SaveTime(save_counter);
        %         end
        %     end
        % end

        % Update Counter
        counter = counter + 1;

        % Store Time, Wellbore Pressure and Injection Rate
        tstore(end+1,1) = t;
        Qstore(end+1,1) = Qt;
        pwstore(end+1,1) = pbot(1); 

    end
    if Control.Manufactured == 1
        m_L2 = ManufacturedSolutionRoutine(hb,pbot,Mesh, Control, t);
        L2store(itern,iter) = m_L2;
    end


    % Plot Injection Rate over Time
    % plt = figure();
    % plot(tstore ./ (60 * 60 * 24 * 365), (Qstore ./ Material.Bg .* (2 * 60 * 60 * 24 * 365) ./ (10^9)))
    % grid on
    % xlabel('Time [yrs]','Interpreter','latex')
    % ylabel('Injection Rate [Mt/yr]','Interpreter','latex')
    % runpath = ['Data ',num2str(iter)];
    % runpath_dir = fullfile(data_dir,runpath,"Injection Rate vs. Time.png");
    % set(gcf, 'PaperPositionMode', 'auto'); % Adjusts the printing to the figure's window size
    % print(plt, runpath_dir, '-dpng', '-r600'); % Saves the figure as a PNG file with 600 DPI resolution
    % runpath_dir = fullfile(data_dir,runpath,"Qsvst.txt");
    % writematrix(num2str([tstore,Qstore]), runpath_dir);












    
    















