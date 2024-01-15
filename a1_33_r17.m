%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mason Shopperly
% AER1410H W Topology Optimization
% Filename: a1_1.m
% Description: Mainline to minimize f under g1 & g2 using Newton's method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockTubeProblem();
function shockTubeProblem()
    clc; close all;

    % 1. Problem and geometry specification
    [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem();
    [L, x0] = specifyGeometry();

    % 2. Mesh Generation
    x_domain = generateMesh(L);

    % 3. Exact solution of governing equations
    [density, mach_number] = solveGoverningEquations(x_domain, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0);

    % 4. Post-processing, assessment, and interpretation of results
    plotResults(x_domain, mach_number, density);
end

function [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem()
    clc; close all;
    % Specify the constants and initial conditions for the shock tube problem
    gamma = 1.4; % Specific heat ratio
    [pL, pR, rhoL, rhoR] = deal(1e5, 1e4, 1, 0.125); % Initial pressures and densities
    [aL, aR] = deal(sqrt(gamma * pL / rhoL), sqrt(gamma * pR / rhoR)); % Speed of sound in left and right sections
    t = 6.1e-3; % Time at which the solution is evaluated
end

function [L, x0] = specifyGeometry()
    L = 10; % Length of the shock tube
    x0 = 5; % Initial location of the diaphragm
end

function x_domain = generateMesh(L)
    % Generate the spatial domain for the shock tube problem
    j = 1000; % Number of points for discretizing the domain
    x_domain = linspace(0, L, j); % Spatial domain
end

function [density, mach_number] = solveGoverningEquations(x_domain, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0)
    % Compute the solution across the domain using the governing equations
    [density, mach_number] = arrayfun(@(x) compute_states(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0), x_domain);
end

function [density, mach_number] = compute_states(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0)
    alpha = (gamma + 1) / (gamma - 1); % Used in several shock and expansion calculations (Equation 3.54)
    
    % The equations are sensitive to the input range for P
    % Ensure the range is physically reasonable, given the problem's context
    P_guess = [0.01 * pR/pL, pL/pR]; % Initial guess for pressure ratio P (Related to Section 3.3)
    
    options = optimset('TolX', 1e-9, 'TolFun', 1e-9);
    % Solve for the pressure ratio using a nonlinear solver
    P = fzero(@(P) P_fun(P, pL, pR, aL, aR, gamma), P_guess, options); % Pressure ratio across the shock (Equation 3.54)
    
    % Calculate the shock speed (C) and contact discontinuity speed (V)
    C = aR * sqrt((gamma + 1)/(2*gamma) * P + (gamma - 1)/(2*gamma)); % Shock speed (Equation 3.58)
    V = 2 * aL / (gamma + 1) * (sqrt(2*gamma/(gamma + 1) * P + (gamma - 1)/(gamma + 1)) - 1); % Speed of contact discontinuity (Equation 3.56)
    
    % Determine the locations of shock and contact discontinuity
    x_shock = x0 + C * t;  % Location of shock to the right of diaphragm (Section 3.3.2)
    x_contact = x0 + V * t; % Location of contact discontinuity follows the shock (Section 3.3.2)
    
    % Calculate the head and tail of the expansion fan
    x_head = x0 - aL * t;  % Head of the expansion fan (Section 3.3.2)
    x_tail = x0 - (aL - V*(gamma + 1)/2) * t;  % Tail of the expansion fan (Section 3.3.2)
    
    % Determine the state based on the x location
    if x < x_head % Region L: Original state to the left of the diaphragm (Equation 3.43)
        density = rhoL;
        mach_number = 0; % Mach number in the original quiescent state
    elseif x >= x_head && x <= x_tail % Region 5: Within the expansion fan (Equation 3.40-3.42)
        [density, mach_number, pressure] = compute_expansion_fan(x, t, x0, aL, pL, rhoL, gamma);
    elseif x > x_tail && x <= x_contact % Region 3: Between the tail of the expansion fan and contact discontinuity (Derived from Section 3.3.2)
        [density, mach_number, pressure] = compute_contact_discontinuity(pL, pR, rhoL, rhoR, gamma, P, V, x_tail, t, x0, aL);
    elseif x > x_contact && x <= x_shock % Region 2: Between contact discontinuity and shock (Derived from Section 3.3.2)
        [density, mach_number, pressure] = compute_shocked_region(pR, rhoR, gamma, P, C, aL);
    else % Region R: Original state to the right of the diaphragm (Equation 3.44)
        density = rhoR;
        mach_number = 0; % Mach number in the original quiescent state
    end
end

function F = P_fun(P, pL, pR, aL, aR, gamma)
    % Solves the implicit equation for pressure P across the shock (Equation 3.54 from textbook)
    alpha = (gamma + 1) / (gamma - 1);  % Definition of alpha (Equation 3.54)
    term1 = sqrt((2 * gamma * (gamma - 1)) / (gamma + 1));  % Part of implicit equation for P (Equation 3.54)
    term2 = (P - 1) / sqrt(1 + alpha * P);  % Part of implicit equation for P (Equation 3.54)
    term3 = (pR / pL)^(1/(2 * gamma)) * (1 - (P / pR)^( (gamma - 1) / (2 * gamma)));  % Part of implicit equation for P (Equation 3.54)
    
    F = term1 * term2 - term3;  % Implicit equation to be solved for P (Equation 3.54)
end

function [density, mach_number, pressure] = compute_expansion_fan(x, t, x0, aL, pL, rhoL, gamma)
    % Computes the state within the expansion fan (Section 3.3.2 and Equations 3.40-3.42)
    u5 = 2 / (gamma + 1) * (aL + (x - x0) / t);  % Velocity within the expansion fan (Derived from Equations 3.40-3.42)
    a5 = aL - (gamma - 1) * u5 / 2;  % Sound speed within the expansion fan (Derived from Equations 3.40-3.42)
    p5 = pL * (a5 / aL)^(2 * gamma / (gamma - 1)); % Pressure within the expansion fan (Derived from Equations 3.40-3.42)
    rho5 = rhoL * (p5 / pL)^(1 / gamma); % Density within the expansion fan (Derived from Equations 3.40-3.42)
    density = rho5;
    velocity = u5;
    mach_number = velocity / a5; % Mach number within the expansion fan (Derived from Equations 3.40-3.42)
    pressure = p5;
end

function [density, mach_number, pressure] = compute_contact_discontinuity(pL, pR, rhoL, rhoR, gamma, P, V, x_tail, t, x0, aL)
    % Computes the state just behind the contact discontinuity (Region 3)
    % This is the state at the tail of the expansion fan, which becomes the state just ahead of the contact discontinuity.
    [~, ~, p_tail] = compute_expansion_fan(x_tail, t, x0, aL, pL, rhoL, gamma); % Pressure at the tail of expansion fan (Used for Region 3)
    density = rhoL * (p_tail / pL)^(1 / gamma);  % Density in Region 3, using isentropic relations (Section 3.3.2)
    velocity = V;  % Velocity is constant across the contact discontinuity (Section 3.3.2)
    mach_number = velocity / aL;  % Mach number just behind the contact discontinuity (Section 3.3.2)
    pressure = p_tail;  % Pressure is continuous across the contact surface (Section 3.3.2)
end

function [density, mach_number, pressure] = compute_shocked_region(pR, rhoR, gamma, P, C, aL)
    % Rankine-Hugoniot conditions to calculate post-shock state
    % Refer to the Rankine-Hugoniot conditions in Section 3.3.2 for the theoretical background
    
    % Calculate post-shock pressure using the pressure ratio P from the implicit equation (3.54)
    pressure = P * pR; % Post-shock pressure (Equation 3.50 in the provided text)
    
    % Calculate post-shock density using equation (3.55) from the provided text
    % This equation is derived from the Rankine-Hugoniot relations
    density = rhoR * ((gamma + 1) * P + (gamma - 1)) / ((gamma + 1) + (gamma - 1) * P); % Post-shock density (Equation 3.55)
    
    % Calculate post-shock velocity using equation (3.58) from the provided text
    % The equation relates the shock speed C to the particle velocity behind the shock
    velocity = (P - 1) * C / (gamma * sqrt((gamma + 1) * P + (gamma - 1))); % Post-shock velocity (Equation 3.58)
    mach_number = velocity / aL; % Mach number using the local speed of sound aL (Equation 3.58 and Mach number definition)
end

function plotResults(x, M, rho)
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