%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
% addpath('2D Array Functions')
% addpath('COCO Continuation/Shapes Point Data/')
% addpath('Shapes Rise Data')
% addpath('Visualize')
% addpath('COCO Continuation/Visualize')
% addpath('COCO Continuation/functions')
% addpath('COCO Continuation')
% addpath('Time Integration')
restoredefaultpath
startup
addpath('..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\General Time Integration Code (MATLAB)\2D Array Functions')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\functions\')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
addpath('..\General Continuation Code (COCO)\Arbitary Shape\')
addpath('Shapes Point Data/')
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 3;

% initialize parameters to sweep over (b, maybe t)
% bpoints = [0.05:0.01:0.20]; %times pi
bpoints = [0.1 0.15]*pi;

%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
shapeNum = 1;
data = init_shape(shapeNum, data);
[data,run_max_E_per_b] = general_COCO(data, bpoints);

% Inside for loop for each b
for b = bpoints
    %% Run time integration for each b
    % Pattern periodic into long chain %Michael
    data = initialize_time_integration(b_val,data);
    

    % Time integration and mitigate edge effects 
    
    % Determine which arch to force/displace %MICHAEL QUESTION
    

    % Time integration and mitigate edge effects
    % Take a look at the initial condition
    plot_system_once(data.A0,data)

    % Prepare for Time Integration
    data.beta = 0.0075;
    T_end = 5000;
    data = determine_coefficient_matrix(data);
    data = determine_starting_vals(data);
    data = determine_modes_to_skip(data);
    data.A0hat = determine_Ahat_from_A(data.A0,data);
    [t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

    %% Look at end of time integration
    plot_system_once(Ahat(end,:)',data)

    data.A0hat = Ahat(end,:)';
    data.A0hat(data.N*data.N_modes+1:end) = 0;

    %% Determine which arch to force/displace %MICHAEL QUESTION

    %data.impose_rotation_at(24) = 1;
    %data.rotation_omega(24) = 2*pi/T_end;
    %data.rotation_mag(24) = 2.5;

    %data.impose_rotation_at(25) = 1;
    %data.rotation_omega(25) = 2*pi/T_end;
    %data.rotation_mag(25) = 2.5;

    % data.impose_rotation_at(2) = 1;
    % data.rotation_omega(2) = 0.0013;
    % data.rotation_mag(2) = 0;
    %
    % data.impose_rotation_at(3) = 1;
    % data.rotation_omega(3) = 0.0013;
    % data.rotation_mag(3) = 0;
    %
    % data.impose_rotation_at(46) = 1;
    % data.rotation_omega(46) = 0.0013;
    % data.rotation_mag(46) = 0;
    %
    % data.impose_rotation_at(47) = 1;
    % data.rotation_omega(47) = 0.0013;
    % data.rotation_mag(47) = 0;
    %
    % data.impose_displacement_at(52) = 0.5;
    % data.displacement_omega(52) = 2*pi/T_end;

    % Remake data with new actuation
    data = determine_coefficient_matrix(data);
    data = determine_starting_vals(data);
    data = determine_modes_to_skip(data);
    %data.A0hat = determine_Ahat_from_A(data.A0,data);

    % We have tiled state with edge effects, now need to compare transition
    % Save data being iterated over
    data_Orig = deepCopyStruct(data);
    % Recover A
    A_Orig = determine_A_from_Ahat(Ahat,data);

    % for loop for different beta values
    %% Run Time Integration for each Beta
    for beta = betavals
        A = A_Orig;
        data = data_Orig;
        %% Prepare for Time Integration
        data.beta = beta;

        %% Run Time Integration
        [t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

        %% Recover A
        A = determine_A_from_Ahat(Ahat,data);

        %% Visualize the Results
        data.frames = 100;
        data.file_name = "b = " + num2str(bval) + " beta = " + num2str(data.beta);
        plot_system_over_time(t,A,data)

        % Determine if a transition occurred and save info (boolean? or distance of wave?)
        % Within each unit cell, calc potential energy, see which ones went from
        % high to low (within say 10% and within 3 unit/super cells)
    end
end