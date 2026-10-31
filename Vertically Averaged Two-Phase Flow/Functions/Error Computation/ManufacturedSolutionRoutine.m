function m_L2 = ManufacturedSolutionRoutine(hb,pbot,Mesh, Control, t)
%   MANUFACTUREDSOLUTION_CHECK Calculates the convergence rates
%   [m_L2, m_e] = ManufacturedSolution_check(hb,Mesh, Control, t)
%   calculates the rates of convergence of the L2 error norm
%
%   ----------------------------------------------------------
%   Acknowledgements: Giovanni Zingaro
%   ----------------------------------------------------------

    % Plot Numerical and Exact Solution
    if Control.plot_toggle == 1
        
        % Presssure Field
        % Compute Exact Solution
        pexact = Control.Manpbot(Mesh.r,t);

        % Create Figure
        figure(1)
        plot(Mesh.r, pbot,'b-.')
        hold on
        plot(Mesh.r, pexact,'r:');
        hold off
        legend('Numerical Solution','Exact Solution','Interpreter','latex')
        grid on
        xlabel('$r$','Interpreter','latex')
        ylabel('$p_{bot}(r)$','Interpreter','latex')

        % Saturation Field
        % Compute Exact Solution
        hbexact = Control.Manhb(Mesh.r,t);

        % Create Figure 
        figure(2)
        plot(Mesh.r, hb,'b-.')
        hold on
        plot(Mesh.r, hbexact,'r:');
        hold off
        legend('Numerical Solution','Exact Solution','Interpreter','latex')
        grid on
        xlabel('$r$','Interpreter','latex')
        ylabel('$h_B(r)$','Interpreter','latex')


    end

        % Compute L2-Error Saturation Field
        if Control.Manufactured == 1
            m_L2 = ManufacturedSolutionCheck(hb,Mesh, Control, t);
        end
end
    