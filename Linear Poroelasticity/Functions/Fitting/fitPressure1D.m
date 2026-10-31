function [p_function] = fitPressure1D(p, disp_fig, r, fittype) 

% Fit Gaussian to Pore Pressure over x and z
    p_function = fit(r,p,fittype);
    p_predict = p_function(r);

    if disp_fig == 1
        close all
        figure(1);
        set(gcf, 'Visible', 'on');

        plot(r,p)
        hold on
        plot(r,p_predict)
        hold off
        legend('Exact','Fitted')
        drawnow
    end
    
end
