%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mason Shopperly
% AER1318H S Topics in Computational Fluid Dynamics
% Filename: a1_31.m
% Description: Mainline to determine the exact solution of the
% subsonic quasi-one-dimensional channel flow problem specified
% in exercise 3.1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a1_31(~, ~)
    clc;
    [gamma, R, T01, p01, S_star] = defineConstants31();
    [j, x] = generateGrid31();
    S = defineAreaFunction31();
    M = solveMachDistribution31(gamma, S, S_star, x);
    [T, p] = calculateDistributions31(M, gamma, T01, p01);
    plotResults31(x, M, p);
    %%
    function [gamma, R, T01, p01, S_star] = defineConstants31()
        gamma = 1.4; % Specific heat ratio
        R = 287; % Universal gas constant [J/(kg*K)]
        T01 = 300; % Total inlet temperature [K]
        p01 = 100e3; % Total inlet pressure [Pa]
        S_star = 0.8; % Critical area for subsonic flow
    end
    %%
    function [j, x] = generateGrid31()
        j = 100;
        Ll = 0;
        Lr = 10;
        x = linspace(Ll, Lr, j);
    end
    %%
    function S = defineAreaFunction31()
        S = @(x) (x >= 0 & x <= 5) .* (1 + 1.5 * (1 - x/5).^2) + ...
                 (x > 5 & x <= 10) .* (1 + 0.5 * (1 - x/5).^2);
    end
    %%
    function [f, df] = areaMachFunction31(M, gamma, Sx, S_star)
        % Calculate the area ratio and its derivative with respect to M
        area_ratio = (1./M) .* ((2/(gamma + 1)) .* (1 + (gamma - 1)/2 .* M.^2)).^((gamma + 1)/(2*(gamma - 1)));
        f = area_ratio - Sx/S_star;
        df = -((gamma + 1)./(2.*(gamma - 1))) .* (1./M.^2) .* ((2./(gamma + 1)) .* (1 + (gamma - 1)/2 .* M.^2)).^((gamma + 1)/(2*(gamma - 1)) - 1) .* ...
             ((2./(gamma + 1)) .* (gamma - 1) .* M);
    end
    %%
    function M = solveMachDistribution31(gamma, S, S_star, x)
        M = zeros(size(x)); % Initialize Mach number array
        M(:) = 0.5; % Initial guess for Mach number (subsonic)
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-9);
        for i = 1:length(x)
            M(i) = fsolve(@(M) areaMachFunction31(M, gamma, S(x(i)), S_star), M(i), options);
        end
    end
    %%
    function [T, p] = calculateDistributions31(M, gamma, T01, p01)
        T = T01 ./ (1 + (gamma - 1)/2 .* M.^2);
        p = p01 * (1 + (gamma - 1)/2 .* M.^2).^(-gamma/(gamma - 1));
    end
    %%
    function plotResults31(x, M, p)
        figure('name', 'Exact solution for the subsonic channel flow problem');
        subplot(2, 1, 1);
        plot(x, p/1000, 'Color', 'k');
        title('Pressure');
        xlabel('x');
        ylabel('p [kPa]');
        yticks(80:2:98);
        ylim([80 98]);
    
        subplot(2, 1, 2);
        plot(x, M, 'Color', 'k');
        title('Mach Number');
        xlabel('x');
        ylabel('M');
        yticks(0.2:.05:.65);
        ylim([0.15 0.65]);
    
        % Set figure position to be centered from top to bottom of the screen
        screen_size = get(0, 'ScreenSize'); % Get the screen size
        fig_width = 400;
        fig_height = 700;
        fig_x = (screen_size(3) - fig_width) / 2; % Horizontally centered
        fig_y = (screen_size(4) - fig_height) / 2; % Vertically centered
        set(gcf, 'Position', [fig_x, fig_y, fig_width, fig_height]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%