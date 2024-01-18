%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mason Shopperly
% AER1318H W Topics in Computational Fluid Dynamics
% Filename: a1_33.m
% Description: Mainline to determine the exact solution of the
% shock-tube problem specified in exercise 3.3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
shockTubeProblem();
function shockTubeProblem()
    clc; close all;

    % 1. Problem & geometry specification
    [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem();
    [x0, x_domain] = specifyGeometry();

    % 2. Exact solution of governing equations & plotting
    [density, mach_number, pressure] = arrayfun(@(x) getState(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0), x_domain);
    plotResults(x_domain, density, mach_number, pressure)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem()
    clc;
    % Initial L & R pressures, densities, specific heat ratio, and solution time as per Exercise 3.3
    [pL, pR, rhoL, rhoR, gamma, t] = deal(1e5, 1e4, 1, 0.125, 1.4, 6.1e-3); 
    % Sound speeds in left and right sections determined from the specified pressures and densities (Equation 3.13)
    [aL, aR] = deal(sqrt(gamma * pL / rhoL), sqrt(gamma * pR / rhoR)); 
end
%%
function [x0, x_domain] = specifyGeometry()
    % Length of the shock tube
    L = 10;
    % Assumed initial location of the diaphragm
    x0 = 5;
    % Number of points for discretizing the domain
    j = 1000; 
    % Spatial domain
    x_domain = linspace(0, L, j);
end
%%
function [density, mach_number, pressure] = getState(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0)
    % Initialization
    [P, V, C, rho3, rho2, p2] = computeShockAndContactSpeeds(pL, pR, aL, aR, rhoL, rhoR, gamma);

    % Define region intersects based on x values
    int_L5 = x0 - aL * t; % Head of the expansion wave
    int_53 = x0 + (V * (gamma + 1) / 2 - aL) * t; % Tail of the expansion fan using correct characteristic speed
    int_32 = x0 + V * t; % Contact surface
    int_2R = x0 + C * t; % Shock wave

    % Determine the state based on the position x
    if x < int_L5
        % Region L
        [density, mach_number, pressure] = computeRegionL(rhoL, pL);
    elseif x >= int_L5 && x < int_53
        % Expansion Fan
        [density, mach_number, pressure] = computeExpansionFan(x, x0, aL, pL, rhoL, gamma, t);
    elseif x >= int_53 && x < int_32
        % Region 3 - Post-Contact
        [density, mach_number, pressure] = computePostContact(p2, rhoL, gamma, pL, V);
    elseif x >= int_32 && x < int_2R
        % Region 2 - Post-Shock
        [density, mach_number, pressure] = computePostShock(P, pR, rhoR, gamma, V);
    else
        % Region R
        [density, mach_number, pressure] = computeRegionR(rhoR, pR);
    end
end
%%
function [P, V, C, rho3, rho2, p2] = computeShockAndContactSpeeds(pL, pR, aL, aR, rhoL, rhoR, gamma)   
    % Modify the initial guesses for pressure ratio P if necessary
    % The interval for P must bracket the root, so we choose P_min and P_max
    % to ensure that they are on opposite sides of the root.
    P_min = pR / pL; % Minimum possible value for P (when the pressure ratio is 1:1)
    P_max = pL / pR; % Maximum possible value for P (when all of the pressure from pL transfers to pR)
    P_guess = [P_min, P_max]; % Updated guess interval
    
    % Function handle for the pressure function
    % This function is derived from the shock jump conditions and the isentropic
    % flow relations for the expansion fan.
    pressureFun = @(P) P_fun(P, pL, pR, aL, aR, gamma);
    
    % Checking if the interval truly brackets a root
    if pressureFun(P_min) * pressureFun(P_max) >= 0
        error('fzero cannot proceed: Function values at the interval endpoints must differ in sign.');
    end
    
    % Solve for the pressure ratio P using a nonlinear solver
    options = optimset('TolX', 1e-9, 'TolFun', 1e-9);
    P = fzero(pressureFun, P_guess, options); 

    % Calculate p2, rho2, V, and C using the obtained pressure ratio P
    % These are derived from the Rankine-Hugoniot relations and the conditions
    % across the contact discontinuity.
    alpha = (gamma + 1) / (gamma - 1); 
    p2 = P * pR; % Pressure behind the shock
    rho2 = rhoR * (1 + alpha * P) / (alpha + P); % Density behind the shock
    p3 = p2; % Pressure is constant across the contact discontinuity
    V = 2 * aL / (gamma - 1) * (1 - (p3 / pL)^((gamma - 1) / (2 * gamma))); % Velocity behind the contact surface
    rho3 = rhoL * (p3 / pL)^(1 / gamma); % Density using isentropic relation
    C = (P - 1) * aR^2 / (gamma * V); % Shock speed (Equation 3.58)
end
%%
function F = P_fun(P, pL, pR, aL, aR, gamma)
    % Implicit equation for pressure P across the shock
    alpha = (gamma + 1) / (gamma - 1);
    term1 = sqrt(2 / (gamma * (gamma - 1))) * (P - 1) / sqrt(1 + alpha * P);
    term2 = (2 / (gamma - 1)) * (aL / aR);
    term3 = (1 - (pR * P / pL)^((gamma - 1) / (2 * gamma)));

    F = term1 - term2 * term3;
end
%%
function [density, mach_number, pressure] = computeRegionL(rhoL, pL)
        density = rhoL;
        mach_number = 0;
        pressure = pL;
end
%%
function [density, mach_number, pressure] = computeExpansionFan(x, x0, aL, pL, rhoL, gamma, t)
    % Compute the local speed of sound within the expansion fan
    u5 = (2 / (gamma + 1)) * ((x - x0) / t + aL);
    a5 = aL - ((gamma - 1) / 2) * u5;

    % Compute the pressure within the expansion fan
    p5 = pL * (a5 / aL)^(2 * gamma / (gamma - 1));

    % Compute the density within the expansion fan
    rho5 = rhoL * (p5 / pL)^(1 / gamma);

    % Compute the Mach number within the expansion fan
    M5 = u5 / a5;

    % Assign outputs to the function
    density = rho5;
    mach_number = M5;
    pressure = p5;
end
%%
function [density, mach_number, pressure] = computePostContact(p2, rhoL, gamma, pL, V)
    % The pressure behind the contact surface is the same as the pressure p2
    % from the post-shock region, and since the flow is isentropic, the density
    % changes according to isentropic relations.
    pressure = p2;
    
    % Calculate the density using isentropic relations, since the entropy to
    % the left of the contact surface is equal to that of the original
    % quiescent left state.
    density = rhoL * (pressure / pL)^(1 / gamma);
    
    % The velocity V is the same across the contact discontinuity and is
    % equal to the speed of the contact surface.
    mach_number = V / (sqrt(gamma * pressure / density)); % Mach number using the speed of sound after the contact surface
end
%%
function [density, mach_number, pressure] = computePostShock(P, pR, rhoR, gamma, V)
    % This function computes the flow properties in region 2, which is 
    % immediately behind the shock wave and is governed by the Rankine-Hugoniot 
    % relations due to the shock.
    
    % Pressure behind the shock is P times the right pressure (pR).
    pressure = P * pR;
    
    % Density behind the shock (rho2) calculated using the shock pressure ratio (P).
    alpha = (gamma + 1) / (gamma - 1);
    density = rhoR * (alpha * P + 1) / (alpha + P);
    
    % Mach number behind the shock is the velocity divided by the speed of sound in region 2.
    mach_number = V / sqrt(gamma * pressure / density);
end
%%
function [density, mach_number, pressure] = computeRegionR(rhoR, pR)
        density = rhoR;
        mach_number = 0;
        pressure = pR;
end
%%
function plotResults(x, rho, M, P)
    close all;
    % This function plots the density and Mach number profiles
    % The figure setup is based on the layout of Fig. 3.3 from the textbook
    
    figure('Name', 'Shock Tube Problem Results');

    subplot(2, 1, 1);
    plot(x, rho, 'Color', 'k');
    title('Density (in Kg/m^3)'); % Title based on Fig. 3.3a in the textbook
    xlabel('x'); % x-axis label as per standard convention
    xlim([0 10]); % x-axis limits as per the domain specified in the textbook's example
    ylabel('\rho', 'Interpreter', 'tex'); % y-axis label using LaTeX formatting as in the textbook
    ylim([0 1.1]); % y-axis limits based on expected range of density values
    yticks(0:.2:1); % y-axis ticks setting for density plot

    subplot(2, 1, 2);
    plot(x, M, 'Color', 'k');
    title('Mach Number'); % Title based on Fig. 3.3b in the textbook
    xlabel('x '); % x-axis label as per standard convention
    xlim([0 10]); % x-axis limits as per the domain specified in the textbook's example
    ylabel('M'); % y-axis label for Mach number
    ylim([0 1]); % y-axis limits based on expected range of Mach number values
    yticks(0:.1:1); % y-axis ticks setting for Mach number plot
    
    % subplot(3, 1, 3);
    % plot(x, P, 'Color', 'k');
    % title('Pressure'); % Title based on Fig. 3.3b in the textbook
    % xlabel('x '); % x-axis label as per standard convention
    % xlim([0 10]); % x-axis limits as per the domain specified in the textbook's example
    % ylabel('P'); % y-axis label for Mach number
    % ylim([1e4 1.1e5]); % y-axis limits based on expected range of Mach number values
    
    % Adjust figure window size and position on screen
    screen_size = get(0, 'ScreenSize');
    fig_width = 400;
    fig_height = 700;
    fig_x = (screen_size(3) - fig_width) / 2;
    fig_y = (screen_size(4) - fig_height) / 2; 
    set(gcf, 'Position', [fig_x, fig_y, fig_width, fig_height]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%