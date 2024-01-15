%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mason Shopperly
% AER1318H W Topics in Computational Fluid Dynamics
% Filename: a1.m
% Description: Mainline to determine the exact solution of the
% quasi-one-dimensional channel flow and the shock-tube problems.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; createGUI();
function flowSolver(caseName)
    global fig;
    global plotS;
    [gamma, p01, S_star, x_shock, j, x] = setupCase(caseName);
    S = defineAreaFunction();
    M = initializeMachNumber(j, x, caseName, x_shock);
    [M, p] = solveFlow(gamma, S, S_star, x, M, p01, caseName, x_shock);
    plotResults(x, M, p, S, caseName, fig, plotS.Value);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supporting functions
function createGUI()
    global plotS;
    global fig;
    f = figure('Position', [100 500 250 250]);
    uicontrol('Style', 'pushbutton', 'String', 'Subsonic',...
              'Position', [20 180 210 20],...
              'Callback', @(~,~) flowSolver('Subsonic'));
    uicontrol('Style', 'pushbutton', 'String', 'Transonic',...
              'Position', [20 150 210 20],...
              'Callback', @(~,~) flowSolver('Transonic'));
    uicontrol('Style', 'pushbutton', 'String', 'Shock-Tube',...
          'Position', [20 120 210 20],...
          'Callback', @(~,~) flowSolver('Shock-Tube'));
    plotS = uicontrol('Style', 'checkbox', 'String', 'Plot S(x)',...
              'Position', [20 90 210 20], 'Value', 0);
    fig = figure('name', 'Exact solution for the flow problems');
    screen_size = get(0, 'ScreenSize');
    fig_width = 400;
    fig_height = 700;
    fig_x = (screen_size(3) - fig_width) / 2;
    fig_y = (screen_size(4) - fig_height) / 2;
    set(fig, 'Position', [fig_x, fig_y, fig_width, fig_height]);
end
function [gamma, p01, S_star, x_shock, j, x] = setupCase(caseName)
    gamma = 1.4;
    R = 287;
    T01 = 300;
    p01 = 100e3;
    if strcmp(caseName, 'Subsonic')
        S_star = 0.8;
        x_shock = NaN;
    elseif strcmp(caseName, 'Transonic')
        S_star = 1;
        x_shock = 7;
    elseif strcmp(caseName, 'Shock-Tube')
        
    end
    j = 1000; % Adjusted number of grid points for resolution
    x = linspace(0, 10, j);
end
function S = defineAreaFunction()
    S = @(x) (x >= 0 & x <= 5) .* (1 + 1.5 * (1 - x/5).^2) + ...
             (x > 5 & x <= 10) .* (1 + 0.5 * (1 - x/5).^2);
end
function [M] = initializeMachNumber(j, x, caseName, x_shock)
    if strcmp(caseName, 'Subsonic')
        M = 0.5 * ones(1, j);
    elseif strcmp(caseName, 'Transonic')
        M = zeros(1, j);
        M(x < 5) = 0.5; 
        M(x >= 5 & x < x_shock) = 2.0; 
        M(x >= x_shock) = 0.5; 
    end
end
function [M, p] = solveFlow(gamma, S, S_star, x, M, p01, caseName, x_shock)
    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
    for i = 1:length(x)
        if strcmp(caseName, 'Subsonic')
            M(i) = fsolve(@(M) areaMachFunction(M, gamma, S(x(i)), S_star), M(i), options);
        elseif strcmp(caseName, 'Transonic')
            if x(i) < x_shock
                M(i) = fsolve(@(M) areaMachFunction(M, gamma, S(x(i)), S_star), M(i), options);
            elseif x(i) == x_shock
                M_post_shock = sqrt((1 + (gamma - 1)/2 * M_pre_shock^2) / (gamma * M_pre_shock^2 - (gamma - 1)/2));
                M(i) = M_post_shock;
            else
                M(i) = fsolve(@(M) areaMachFunction(M, gamma, S(x(i)), S_star), M(i), options);
            end
        end
        p(i) = p01 * (1 + (gamma - 1)/2 * M(i)^2)^(-gamma/(gamma - 1));
    end
end
function F = areaMachFunction(M, gamma, Sx, S_star)
    % Equation 3.45
    F = (Sx/S_star) - (1/M) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * M^2))^((gamma + 1)/(2*(gamma - 1)));
end
function plotResults(x, M, p, S, caseName, fig, plotS)
    figure(fig);
    clf(fig, 'reset'); 
    subplot(2, 1, 1);
    hold on;
    plot(x, p/1000, 'k', 'DisplayName', 'Pressure');
    ylabel('p [kPa]');
    if plotS
        yyaxis right;
        plot(x, arrayfun(S, x), 'r:', 'DisplayName', 'S(x)', 'LineWidth', 1.5);
        ylabel('S(x)');
        ylim([0.5 3]);
        legend('Location', 'northeast');

    end
    title(['Pressure - ' caseName]);
    xlabel('x');
    hold off;
    subplot(2, 1, 2);
    hold on;
    plot(x, M, 'k', 'DisplayName', 'Mach Number');
    ylabel('M');
    if plotS
        yyaxis right;
        plot(x, arrayfun(S, x), 'r:', 'DisplayName', 'S(x)', 'LineWidth', 1.5);
        ylabel('S(x)');
        ylim([0.5 3]);
        legend('Location', 'northeast');        
    end
    title(['Mach Number - ' caseName]);
    xlabel('x');
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%