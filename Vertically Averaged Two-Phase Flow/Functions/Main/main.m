%% Axisymmetric Vertically Averaged Multi-Phase Flow Code -> Main File
% Acknowledgements: Giovanni Zingaro

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
    h.hcn1 = (Mesh.H)*ones(Mesh.nn,1)-h.hbn1;
    h.hc = h.hcn1;
    h.hb = h.hbn1;

%% Storage Arrays
    tstore = [];
    Qstore = [];
    pwstore = [];

%% Initialize Time and Spatial Parameters
    t = Control.StartTime; % Current Time
    tend = Control.EndTime; % End Time
    theta = Control.theta; % Numerical Integration Parameter
    dr = Mesh.L/Mesh.ner; % Spatial Discretization

    % Store Time, Wellbore Pressure and Injection Rate
    tstore(end+1,1) = 0;
    Qstore(end+1,1) = 0;

    % Initialize Initial Injection Rate
    Qti = Qmean_i + Qstd_i * randn(1);
    Qt = max(Qmin,Qti)* Material.Bg.*(1e9)./(2 * 60 * 60 * 24 * 365);

    tstore(end+1,1) = 0;
    Qstore(end+1,1) = Qt;

    % Update Boundary Condition for Flux of CO2
    if Control.Manufactured == 0
        BC.qc = @(r,t) Q0;
    end

%% Assemble and Solve Pressure Equations for Initial Condition
    
    % Linear Diffusion/Stiffness Matrix for Pressure Equation
    Kp = getKp(Mesh, Quad, Material, h); % Diffusion/Stiffness Matrix

    % Injection Flux Vector for Pressure Equation
    Fp2 = getFp2(Mesh, Quad, BC, t);

    % Buonyancy Vector for Pressure Equation
    Fp3 = getFp3(Mesh, Quad, h.hc, Material);

    % Sum the External Flux Vector
    Fp = Fp2 - Fp3;

    % Solve the System of Equations for the Pressure Equation
    [pbot,~] = LinearSolver1(Kp,Fp,BC.fix_pres_value,BC.free_pres,BC.fixed_pres,Control.parallel);

    % Save Wellbore Pressure
    pwstore(end+1,1) = pbot(1); 

%% Declare Timestep for First Step
    dtn1 = Control.TimeStep;
    counter = 1; % Timestep Counter 

%% Save Intial Conditions

    % Save Counter 
    save_counter = 1;
    counter_end = 0;

    % Save Initial Conditions in Textfiles
    % if Control.Save == 1
    %     writeData2txt(pbot, h.hb, Mesh, Material, data_dir, iter, Qt, t, counter, BC, counter_end);
    % end

    if (Control.Save) == 1

        % Next Saving Time
        save_counter = save_counter + 1;
        tsave = Control.SaveTime(save_counter);

    end

%% Declare Mass Conservation Trackers
Qex = 0;

%% Timestep Loop
    while t<=tend

        %% Assemble and Solve Pressure Equations

        % Update Boundary Condition for Flux of CO2
        if Control.Manufactured == 0
            BC.qc = @(r,t) Qt;
        end

        % Linear Diffusion/Stiffness Matrix for Pressure Equation
        Kp = getKp(Mesh, Quad, Material, h); % Diffusion/Stiffness Matrix

        % Boundary Flux Vector for Pressure Equation
        Fp1 = getFp1(Mesh, t, BC); 

        % Injection Flux Vector for Pressure Equation
        Fp2 = getFp2(Mesh, Quad, BC, t);

        % Buonyancy Vector for Pressure Equation
        Fp3 = getFp3(Mesh, Quad, h.hc, Material);

        % Sum the External Flux Vector
        Fp = Fp1 + Fp2 - Fp3;

        % Solve the System of Equations for the Pressure Equation
        [pbot,~] = LinearSolver1(Kp,Fp,BC.fix_pres_value,BC.free_pres,BC.fixed_pres,Control.parallel);

        % Compute Minimum Timestep
        dt = ComputeTimeStep(Mesh, Quad, Material, pbot, Control, BC);
        dt = min(dt, tsave - t);

        %% Modify Time and Injection Rates
        if t == tend
            break
        elseif t+dt > tend
            dt = tend-t;
            t = t + dt;

        elseif (t+dt) >= tsave
            
                % Modify Timestep
                dt = tsave-t;
                
                % Update Time
                t = t + dt;

                % Save Time and Injection Rate
                tstore(end+1,1) = t;
                Qstore(end+1,1) = Qt;

        else

            % Update Time 
            t = t + dt;
        end

        % Declare Adpative Timestep
        Control.TimeStep = dt;

        %% Solve Saturation Equation
        % Assemble Fully Stabilized Advection Reaction Matrix (K)
        Ks = getKsu(Mesh, Quad, Material, pbot, Control);

        % Declare Pre-Computed Storage Matrix
        C = getCsu(Mesh, Quad); % Fully Stabilized Storage Matrix

        % Declare Flux Vector for Saturation Equation
        Fg = zeros(Mesh.nn,1);

        %% Solve Saturation Equation

        % Solve the Transient Problem
        Keff = C+dt*theta*Ks;
        Feff = (C-dt*(1-theta)*Ks)*(h.hbn1)+dt*(theta*Fg+(1-theta)*Fg);

        % Solve the System of Equations for the Saturation Equation
        [hb,fE] = LinearSolver1(Keff,Feff,BC.fixed_sat_value,BC.free_sat,BC.fixed_sat,Control.parallel);

        % Constraint Limiter 
        negativeIndex = find(hb < 0);

        % Check if a negative value was found
        if ~isempty(negativeIndex)

            % Set all values from the first negative index to the end of the array to 0
            hb(negativeIndex) = 0;

        end

        % Constraint Limiter 
        positiveIndex = find(hb > Mesh.H);

        % Check if a negative value was found
        if ~isempty(positiveIndex)

            % Set all values from the first negative index to the end of the array to 0
            hb(positiveIndex) = Mesh.H;

        end

        hb(2000:end) = Mesh.H;

        % Compute Saturation of CO2 
        h.hc = BC.H - hb;

        % Total CO2 Injected 
        Qex = Qex + Qt*dt;
         
        % Update hb and hc
        h.hbn1 =  BC.H - h.hc;
        h.hcn1 = h.hc;

        % Save to Textfile
        if Control.Save == 1
            if t == tsave

                % Scale CO2 Saturation to conserve mass 
                % Check Mass Conservation
                [~, mass_error] = MassError(Mesh, h.hc, Qex, Material);
                hc = h.hc; 

                % % Modify Increment of CO2 Saturation to Conserve Mass
                if mass_scale == 1 && mass_error > (Control.tol)
                    while mass_error > (Control.tol)
                        [dhc_scale, mass_error] = MassError(Mesh, round(hc,8), Qex, Material); % Compute Mass Error and Scale
                        hc = dhc_scale.*hc; % Update CO2 Saturation
                    end
                end

                % Update Brine Saturation for Saving
                hb = BC.H - hc;

                % Plot Saturation at End
                if t == tend
                    counter_end = 1;
                end

                % Save data to text file
                writeData2txt(pbot, hb, Mesh, Material, data_dir, iter, Qt, t, counter, BC, counter_end);

                % Modify Flow Rate
                Qti = Qmean_i + Qstd_i * randn(1);
                Qt = max(Qmin,Qti)* Material.Bg.*(1e9)./(2 * 60 * 60 * 24 * 365);

                % Store Time and Injection Rate
                tstore(end+1,1) = t;
                Qstore(end+1,1) = Qt;
                
                if t == Control.SaveTime(end)
                else
                    save_counter = save_counter + 1;
                    tsave = Control.SaveTime(save_counter);
                end
            end
        end

        % Update Counter
        counter = counter + 1;

        % Store Time, Wellbore Pressure and Injection Rate
        pwstore(end+1,1) = pbot(1); 

    end
    if Control.Manufactured == 1
        m_L2 = ManufacturedSolutionRoutine(hb,pbot,Mesh, Control, t);
        L2store(itern,iter) = m_L2;
    end

    % Save data to text file
    %counter_end = 1;
    %writeData2txt(pbot, hc, Mesh, Material, data_dir, iter, Qt, t, counter, BC, counter_end);

    % Plot Injection Rate over Time
    plt = figure(Visible="off");
    plot(tstore ./ (60 * 60 * 24 * 30), (Qstore ./ Material.Bg .* (2 * 60 * 60 * 24 * 365) ./ (10^9)))
    grid on
    xlabel('Time [Months]','Interpreter','latex')
    ylabel('Injection Rate [Mt/yr]','Interpreter','latex')
    runpath = ['Data ',num2str(iter)];
    runpath_dir = fullfile(data_dir,runpath,"Injection Rate vs. Time.png");
    set(gcf, 'PaperPositionMode', 'auto'); % Adjusts the printing to the figure's window size
    print(plt, runpath_dir, '-dpng', '-r600'); % Saves the figure as a PNG file with 600 DPI resolution
    runpath_dir = fullfile(data_dir,runpath,"Qsvst.txt");
    writematrix(num2str([tstore,Qstore]), runpath_dir);












    
    















