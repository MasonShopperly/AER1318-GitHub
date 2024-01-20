%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mason Shopperly
% AER1318H S Topics in Computational Fluid Dynamics
% Filename: a1_32.m
% Description: Mainline to determine the exact solution of the
% transonic quasi-one-dimensional channel flow problem specified
% in exercise 3.2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a1_32(~, ~)
    clc;
    [gamma, R, T01, p01, S_star, x_shock, j, x, rho01, a01, pL_star_aL_star] = setupSimulation32();
    S = defineAreaFunction32();
    [M] = initializeMachNumber32(j, x);
    [M, p] = solveTransonicFlow32(gamma, S, S_star, x, M, p01, R, T01, x_shock);
    plotResults32(x, M, p);
    %%
    function [gamma, R, T01, p01, S_star, x_shock, j, x, rho01, a01, pL_star_aL_star] = setupSimulation32()
        gamma = 1.4; % Specific heat ratio
        R = 287; % Universal gas constant [J/(kg*K)]
        T01 = 300; % Total inlet temperature [K]
        p01 = 100e3; % Total inlet pressure [Pa]
        S_star = 1; % Critical area
        x_shock = 7; % Shock location
        j = 1000; % Number of grid points
        x = linspace(0, 10, j); % Grid space

        rho01 = p01 / (R * T01); 
        a01 = sqrt(gamma * p01 / rho01);
        pL_star_aL_star = p01 * a01 * (2 / (gamma + 1))^((gamma  + 1) / (2 *( gamma + 1)));
    end
    %%
    function S = defineAreaFunction32()
        S = @(x) (x >= 0 & x <= 5) .* (1 + 1.5 * (1 - x/5).^2) + ...
                 (x > 5 & x <= 10) .* (1 + 0.5 * (1 - x/5).^2);
    end
    %%
    function [M] = initializeMachNumber32(j, x)
        M = zeros(1, j);
        for i = 1:j
            if x(i) < 5
                M(i) = 0.5; 
            elseif x(i) >= 5 && x(i) < 7
                M(i) = 1.2;
            else
                M(i) = 0.6;
            end
        end
    end
    %%
    function [M, p] = solveTransonicFlow32(gamma, S, S_star, x, M, p01, R, T01, x_shock)
        % Solve for Mach number distribution and calculate pressure
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-9);
        shock_index = find(x >= x_shock, 1); % Find the index where the shock occurs
        
        % Initialize variables
        p = zeros(1, length(x)); % Pressure
        S_star_R = S_star; % Start with pre-shock area ratio
        p0_R = p01; % Start with pre-shock total pressure
    
        for i = 1:length(x)
            if i < shock_index
                % Solve for Mach number using fsolve for pre-shock conditions
                M(i) = fsolve(@(M) areaMachRelation32(M, gamma, S(x(i)), S_star), M(i), options);
                p(i) = p01 * (1 + (gamma - 1) / 2 * M(i)^2)^(-gamma/(gamma - 1));
            elseif i == shock_index
                % Calculate post-shock conditions using Rankine-Hugoniot relations at the shock
                M_L = M(i-1); % Mach number just before the shock
                p_L = p(i-1); % Pressure just before the shock
                
                % Calculate downstream Mach number and pressure using Rankine-Hugoniot relations
                M_R = sqrt((2 + (gamma - 1) * M_L^2) / (2 * gamma * M_L^2 - (gamma - 1)));
                p_R = (2 * gamma * M_L^2 - (gamma - 1)) / (gamma + 1) * p_L;
                
                % Recalculate post-shock total pressure and total density
                p0_R = p_R * (1 + (gamma - 1) / 2 * M_R^2)^(gamma / (gamma - 1));
                rho0_R = p0_R / (R * T01);
                
                % Calculate the post-shock area ratio S_star_R
                a0_L = sqrt(gamma * R * T01); % Assuming T0L = T01, since T0 is constant across a shock
                a0_R = sqrt(gamma * p0_R / rho0_R);
                S_star_R = S_star * (p01 * a0_L) / (p0_R * a0_R);
                
                % Update Mach number and pressure at the shock
                M(i) = M_R;
                p(i) = p_R;
            else
                % For areas downstream of the shock, use post-shock total pressure and area ratio
                M(i) = fsolve(@(M) areaMachRelation32(M, gamma, S(x(i)), S_star_R), M(i-1), options);
                p(i) = p0_R * (1 + (gamma - 1) / 2 * M(i)^2)^(-gamma/(gamma - 1));
            end
        end
    end
    %%
    function F = areaMachRelation32(M, gamma, Sx, S_star)
        % Area-Mach function based on Equation (3.45)
        F = (Sx / S_star) - (1 / M) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M^2))^((gamma + 1) / (2 * (gamma - 1)));
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
