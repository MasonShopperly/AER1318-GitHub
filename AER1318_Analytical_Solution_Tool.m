%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mason Shopperly
% AER1318H W Topics in Computational Fluid Dynamics
% Filename: AER1318_Analytical_Solution_Tool.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AER1318_Analytical_Solution_Tool()
    % Create a figure for the GUI
    fig = figure('Name', 'AER1318 Analytical Solution Tool', ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 550, 350, 200]);

    % Button for Exercise 3.1
    uicontrol('Style', 'pushbutton', 'String', 'Exercise 3.1', ...
              'Position', [50, 150, 250, 30], ...
              'Callback', @a1_31);

    % Button for Exercise 3.2
    uicontrol('Style', 'pushbutton', 'String', 'Exercise 3.2', ...
              'Position', [50,110, 250, 30], ...
              'Callback', @a1_32);

    % Button for Exercise 3.3
    uicontrol('Style', 'pushbutton', 'String', 'Exercise 3.3', ...
              'Position', [50, 70, 250, 30], ...
              'Callback', @a1_33);

    % Button to Reset Plots
    uicontrol('Style', 'pushbutton', 'String', 'Reset Plots', ...
              'Position', [50, 30, 250, 30], ...
              'Callback', @resetPlots);

    function resetPlots(~, ~)
        figs = findall(groot, 'Type', 'figure');
        gui_fig = gcf;
        close(setdiff(figs, gui_fig));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%