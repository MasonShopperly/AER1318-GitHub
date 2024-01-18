%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mason Shopperly
% AER1318H W Topics in Computational Fluid Dynamics
% Filename: a1_32.m
% Description: Mainline to determine the exact solution of the
% transonic quasi-one-dimensional channel flow problem specified
% in exercise 3.2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a1_32(~, ~)
    clc;
    [gamma, R, T01, p01, S_star, x_shock, j, x] = setupSimulation32();
    S = defineAreaFunction32();
    [M, M_post_shock] = initializeMachNumber32(j, x, gamma);
    [M, p] = solveTransonicFlow32(gamma, S, S_star, x, M, p01);
    plotResults32(x, M, p);
    %%
    function [gamma, R, T01, p01, S_star, x_shock, j, x] = setupSimulation32()
        gamma = 1.4; % Specific heat ratio
        R = 287; % Universal gas constant [J/(kg*K)]
        T01 = 300; % Total inlet temperature [K]
        p01 = 100e3; % Total inlet pressure [Pa]
        S_star = 1; % Critical area
        x_shock = 7; % Shock location
        j = 1000; % Number of grid points
        x = linspace(0, 10, j); % Grid space
    end
    %%
    function S = defineAreaFunction32()
        S = @(x) (x >= 0 & x <= 5) .* (1 + 1.5 * (1 - x/5).^2) + ...
                 (x > 5 & x <= 10) .* (1 + 0.5 * (1 - x/5).^2);
    end
    %%
    function [M, M_post_shock] = initializeMachNumber32(j, x, gamma)
        M = zeros(1, j);
        M_post_shock = NaN; 
        for i = 1:j
            if x(i) < 5
                % Use a subsonic guess before the shock
                M(i) = 0.5; 
            elseif x(i) >= 5 && x(i) < 7
                % Use a supersonic guess before the shock
                M(i) = 1.2;
            else
                % If shock has been passed, use Rankine-Hugoniot relations for post-shock Mach number
                if isnan(M_post_shock)
                    M_post_shock = calculatePostShockMach32(M(i-1), gamma);
                end
                % Use a subsonic guess after the shock
                M(i) = 0.6;
            end
        end
    end
    %%
    function [M, p] = solveTransonicFlow32(gamma, S, S_star, x, M, p01)
        % Solve for Mach number distribution and calculate pressure
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6);
        for i = 1:length(x)
            % Iteratively solve for Mach number using fsolve
            M(i) = fsolve(@(M) areaMachRelation32(M, gamma, S(x(i)), S_star), M(i), options);
            % Calculate pressure using isentropic relations
            p(i) = p01 * (1 + 0.5 * (gamma - 1) * M(i)^2)^(-gamma/(gamma-1));
        end
    end
    %%
    function F = areaMachRelation32(M, gamma, Sx, S_star)
        % Correct implementation of Equation 3.45
        F = (Sx/S_star) - (1/M) * ((2/(gamma + 1)) * (1 + (gamma - 1)/2 * M^2))^((gamma + 1)/(2*(gamma - 1)));
    end
    %%
    function M_post_shock = calculatePostShockMach32(M_pre_shock, gamma)
        % Correct implementation of Rankine-Hugoniot relations (Equation 3.49)
        M_post_shock = sqrt((1 + (gamma - 1)/2 * M_pre_shock^2) / (gamma * M_pre_shock^2 - (gamma - 1)/2));
    end
    %%
    function plotResults32(x, M, p)
        figure;
        subplot(2, 1, 1);
        plot(x, p, 'Color', 'k');
        title('Pressure');
        xlabel('x');
        ylabel('p [Pa]');
        subplot(2, 1, 2);
        plot(x, M, 'Color', 'k');
        title('Mach Number');
        xlabel('x ');
        ylabel('M');
        ylim([.2 1.6]);
        screen_size = get(0, 'ScreenSize');
        fig_width = 400;
        fig_height = 700;
        fig_x = (screen_size(3) - fig_width) / 2;
        fig_y = (screen_size(4) - fig_height) / 2; 
        set(gcf, 'Position', [fig_x, fig_y, fig_width, fig_height]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
