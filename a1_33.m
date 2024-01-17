%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mason Shopperly
% AER1318H W Topics in Computational Fluid Dynamics
% Filename: a1_33.m
% Description: Mainline to determine the exact solution of the
% shock-tube problem specified in exercise 3.3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockTubeProblem();
function shockTubeProblem()
    clc; close all;

    % 1. Problem and geometry specification
    [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem();
    [x0, x_domain] = specifyGeometry();

    % 2. Exact solution of governing equations
    [density, mach_number] = arrayfun(@(x) getState(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0), x_domain);

    % 3. Post-processing, assessment, and interpretation of results
    plotResults(x_domain, mach_number, density);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gamma, pL, pR, rhoL, rhoR, aL, aR, t] = specifyProblem()
    clc;
    % Initial L & R pressures, densities, specific heat ratio, and solution time as per Exercise 3.3
    [pL, pR, rhoL, rhoR, gamma, t] = deal(1e5, 1e4, 1, 0.125, 1.4, 6.1e-3); 
    % Sound speeds in left and right sections determined from the specified pressures and densities (Equation 3.13)
    [aL, aR] = deal(sqrt(gamma * pL / rhoL), sqrt(gamma * pR / rhoR)); 
end

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

function [density, mach_number] = getState(x, t, pL, pR, rhoL, rhoR, aL, aR, gamma, x0)
    % (Equation 3.54)
    alpha = (gamma + 1) / (gamma - 1); 
    
    % Initial guess for pressure ratio P (Related to Section 3.3)
    P_guess = [0.01 * pR/pL, pL/pR]; 
    
    options = optimset('TolX', 1e-9, 'TolFun', 1e-9);
    % Solve for the pressure ratio using a nonlinear solver
    P = fzero(@(P) P_fun(P, pL, pR, aL, aR, gamma), P_guess, options); % Pressure ratio across the shock (Equation 3.54)
    %P = newtonsMethod(pL, pR, aL, aR, gamma, P_guess);
    
    % Calculate the shock speed (C) and contact discontinuity speed (V)
    C = aR * sqrt((gamma + 1)/(2*gamma) * P + (gamma - 1)/(2*gamma)); % Shock speed (Equation 3.58)
    V = 2 * aL / (gamma + 1) * (sqrt(2*gamma/(gamma + 1) * P + (gamma - 1)/(gamma + 1)) - 1); % Speed of contact discontinuity (Equation 3.56)
    
    % Determine the locations of shock and contact discontinuity
    x_shock = x0 + C * t;  % Location of shock to the right of diaphragm (Section 3.3.2)
    x_contact = x0 + V * t; % Location of contact discontinuity follows the shock (Section 3.3.2)
    
    % Calculate the head and tail of the expansion fan
    % Calculate the head of the expansion fan
    x_head = x0 - aL * t;  % Head of the expansion fan moves to the left

    % Calculate the tail of the expansion fan using the correct formula
    u5 = aL - ((gamma + 1) / 2) * V;
    x_tail = x0 - u5 * t;  % Tail of the expansion fan moves to the left
    
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
    % Debugging output
    disp(['x: ', num2str(x), ' t: ', num2str(t)]);
    disp(['x_shock: ', num2str(x_shock), ' x_contact: ', num2str(x_contact)]);
    disp(['x_head: ', num2str(x_head), ' x_tail: ', num2str(x_tail)]);
end

function P = newtonsMethod(pL, pR, aL, aR, gamma, P_guess)
    % Constants for Newton's method
    maxIter = 10000;      % Maximum number of iterations
    tol = 1e-6;         % Tolerance for convergence
    P = P_guess;        % Initial guess for P

    for iter = 1:maxIter
        F = P_fun(P, pL, pR, aL, aR, gamma);   % Evaluate function at current P
        dF = dP_fun(P, pL, pR, aL, aR, gamma); % Evaluate derivative at current P

        % Check if derivative is too close to zero
        if abs(dF) < eps
            error('Derivative too small, stopping iterations');
        end
        
        % Newton-Raphson update
        deltaP = F / dF;
        
        % Dynamic step size adjustment if the Newton step is too large
        while abs(deltaP) > 0.5 * abs(P)
            deltaP = deltaP / 2;
        end
        
        P_new = P - deltaP;

        % Check for convergence
        if abs(F) < tol
            P = P_new;
            break;
        end
        
        P = P_new; % Update P for the next iteration
    end

    % If Newton's method did not converge, throw an error
    if iter == maxIter
        error('Newton''s method did not converge within the maximum number of iterations');
    end
end

function F = P_fun(P, pL, pR, aL, aR, gamma)
    % Solves the implicit equation for pressure P across the shock (Equation 3.54 from textbook)
    alpha = (gamma + 1) / (gamma - 1);  % Definition of alpha (Equation 3.54)
    term1 = sqrt((2 * gamma * (gamma - 1)) / (gamma + 1));  % Part of implicit equation for P (Equation 3.54)
    term2 = (P - 1) / sqrt(1 + alpha * P);  % Part of implicit equation for P (Equation 3.54)
    
    % Ensure P/pR is positive to avoid complex numbers
    if P / pR <= 0
        error('P/pR must be positive to avoid complex results in P_fun.');
    end
    term3 = (pR / pL).^(1/(2 * gamma)) * (1 - (P / pR).^((gamma - 1) / (2 * gamma)));  % Part of implicit equation for P (Equation 3.54)
    
    F = term1 * term2 - term3;  % Implicit equation to be solved for P (Equation 3.54)
end

function dF = dP_fun(P, pL, pR, aL, aR, gamma)
    % Derivative of the function P_fun with respect to P
    alpha = (gamma + 1) / (gamma - 1);
    term1 = sqrt((2 * gamma * (gamma - 1)) / (gamma + 1));
    term2 = 1 ./ sqrt(1 + alpha * P);
    term3 = (alpha / 2) * (P - 1) / ((1 + alpha * P).^(3 / 2));

    % Ensure P/pR is positive to avoid complex numbers
    if P / pR <= 0
        error('P/pR must be positive to avoid complex results in dP_fun.');
    end
    term4 = (pR / pL)^(1 / (2 * gamma));
    term5 = ((gamma - 1) / (2 * gamma)) * (P / pR).^((-1 + gamma) / (2 * gamma));
    
    % Calculate derivative using the power rule
    dP_dTerm3 = -(pR / pL)^(1 / (2 * gamma)) * ((gamma - 1) / (2 * gamma)) * ...
                (P / pR).^((-1 - gamma) / (2 * gamma)) / pR;
    
    dF = term1 * (term2 - term3) - term4 * dP_dTerm3;
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
    % Correct computation of pressure in Region 3
    pressure = P * pR; % Pressure in Region 3 should match that of Region 2
    
    % Correct computation of density in Region 3 using isentropic relation
    density = rhoL * (pressure / pL)^(1 / gamma);
    
    % Velocity in Region 3 is the same as that in Region 2, which is V
    velocity = V;  % Constant across the contact discontinuity
    
    % Recalculate the local speed of sound in Region 3 using the updated pressure and density
    a3 = sqrt(gamma * pressure / density);
    
    % Calculation of the Mach number using the velocity and local speed of sound in Region 3
    mach_number = velocity / a3;
    
    % Debugging output
    disp(['Region 3 - Contact Discontinuity: Pressure = ', num2str(pressure), ...
          ', Density = ', num2str(density), ', Velocity = ', num2str(velocity), ...
          ', Mach Number = ', num2str(mach_number), ...
          ', Speed of Sound = ', num2str(a3)]);
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