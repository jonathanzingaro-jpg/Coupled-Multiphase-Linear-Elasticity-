function writeData2txt(pbot, hb, Mesh, Material, data_dir, iter, Qt, t, counter, BC, counter_end)
%writeData2txt Exports results to data as a text file located in the
%   data_dir
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
%   pbot:       Column vector of pore pressures (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   hb:         Column vector of brine saturations (size ndof x 1)
%   Material:   Structure containing material parameters
%   BC:         Structure containing node sets
%   data_dir:   Directory where the data files are to be saved
%   iter:       Iteration number
%   Qt:         Total flow rate
%   t:          Time in seconds
%   counter:    Counter for determining whether to export additional data
%
%   Acknowledgement: Jonathan Zingaro

    % Assuming injection_rate and time variables are defined
    injection_rate = Qt;  % Injection rate in m^3/day

    % Header information
    header = ["Qs [m^3/s]: " + num2str(injection_rate, '%.16E'), ...
              "Time [days]: " + num2str(t/24/60/60, '%.12f'), ...
              "k [mD]: " + num2str(Material.kx/9.869233e-16, '%.12f'), ...
              "phi [-]: " + num2str(Material.phi, '%.12f'), ...
              "--------------------------------------"]';

    % Convert displacement values to strings with specified format
    pres_str = compose("%.16E", pbot);
    
    % Combine header and formatted displacement strings
    pres_data = [header; pres_str];
    
    % Export Pore Pressure, Saturation, Time, Injection Rate, Reservoir Properties
    runpath = ['Data ',num2str(iter)];
    runpath_dir = fullfile(data_dir, runpath, ['t = ',num2str(round(t/24/60/60,1)), ' days']);
    mkdir(runpath_dir);
    writematrix(pres_data, fullfile(runpath_dir, 'pres.txt'));

    % Similarly for other files
    % You need to convert the numeric data to string to concatenate with the header
    if counter == 1
        nodal_coordinates_data = [header; string(Mesh.r)];
        writematrix(nodal_coordinates_data, fullfile(runpath_dir, 'Nodal_Coordinates.txt'),'Delimiter','comma');
    end

    saturation_brine_data = [header; string(compose("%.16E",hb))];
    writematrix(saturation_brine_data, fullfile(runpath_dir, 'sat_brine.txt'));

    if counter == 1
        % Export Input Functions
        writetable(struct2table(Material), fullfile(runpath_dir,'Flow Properties.csv'))
    end

    %% Plot Pressure and Saturation over Time
    if counter_end == 1
        set(0, 'DefaultFigureVisible', 'off');
        plt = figure;
        subplot(2,1,1)
        plot(Mesh.r,pbot./1e6)
        xlim([0,max(Mesh.r)])
        xlabel('$r [m]$', 'Interpreter', 'latex'); % LaTeX for x-label
        ylabel('$p_{bot} [MPa]$', 'Interpreter', 'latex'); % LaTeX for y-label
        title(['Pressure Build-up at Bottom ($t = ', num2str(t/24/60/60, '%0.2f'), '$ days)'], 'Interpreter', 'latex')
        ylim([0,ceil(max(pbot)/1e6)+1])
        grid on

        subplot(2,1,2)
        plot(Mesh.r, hb)
        xlabel('$r [m]$', 'Interpreter', 'latex'); % LaTeX for x-label
        ylabel('$h_{B} [m]$', 'Interpreter', 'latex'); % LaTeX for y-label
        xlim([0,max(Mesh.r)])
        ylim([0,max(hb)*1.1])
        grid on

        % Add a horizontal line at hb = BC.H
        line(get(gca, 'xlim'), [BC.H BC.H], 'Color', 'black', 'LineStyle', '--');

        % Format the title with LaTeX
        title(['Brine Saturation ($t = ', num2str(t/24/60/60, '%0.2f'), '$ days)'], 'Interpreter', 'latex')

        % Save Pressure and Saturation Distribution at End Time
        filename = fullfile(data_dir, runpath, ['t = ',num2str(round(t/24/60/60,1)), ' days'], ['t = ', num2str(round(t/24/60/60, 1)), ' days', '.png']); % This creates a full file path with .png extension
        set(gcf, 'PaperPositionMode', 'auto'); % Adjusts the printing to the figure's window size
        print(plt, filename, '-dpng', '-r150'); % Saves the figure as a PNG file with 600 DPI resolution
    end
end
