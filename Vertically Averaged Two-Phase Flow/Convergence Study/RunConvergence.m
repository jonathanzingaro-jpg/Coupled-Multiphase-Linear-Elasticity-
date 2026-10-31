%% Initialization and Configuration
clearvars -global
clear, clc, close all
format compact
tic;

%% Define Number of Element and Stability Constants
n = [10 50 100 500 1000];  % Number of elements
dx_store = 1./n;
element_type = {'L2'};
legendLabels = {};  % Initialize legend labels

%% Simulation Loop
global N L2store iter etype itern

L2store = zeros(length(element_type), length(n));
for j = 1:length(element_type)
    itern = j;
    etype = element_type{j};
    for iter = 1:length(n)
        fprintf('\n');
        disp(['Simulation ', num2str(iter)]);
        N = n(iter);  % Number of element in r direction

        run('RunMan.m');  % SU method script
        disp('Done');
    end

    % Compute Slope
    pL2 = polyfit(log(dx_store), log(L2store(itern,:)), 1);
    slopeStr = sprintf('m = %.3f', pL2(1));
    disp(slopeStr);
    legendLabels{end+1} = sprintf('%s: m = %.3f', etype, pL2(1));

    % L2-Error Plot
    loglog(dx_store, L2store(itern,:), '-o','DisplayName', sprintf('%s: %s', etype, slopeStr));
    xlabel('Element Size', 'Interpreter', 'latex');
    ylabel('$L_{2}$-Error', 'Interpreter', 'latex');
    grid on;
    hold on;

end

%% Final Plot Adjustments
legend(legendLabels, 'Location', 'best');
hold off;
drawnow;
