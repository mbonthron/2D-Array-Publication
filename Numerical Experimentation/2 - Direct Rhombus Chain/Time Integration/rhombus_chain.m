%% Clear Everything so there are no stragglers
clear; clc; close all

betavals = .00002;
bvals = [.05:.05:.30];

%% Add the Paths to the Required Functions
restoredefaultpath
startup
addpath('..\..\..\General Time Integration Code (MATLAB)\2D Array Functions\')
addpath('2D Array Functions\')
addpath('Shapes Point Data')
addpath('Shapes Rise Data')
addpath('..\..\..\General Time Integration Code (MATLAB)\Visualize\')

for bval = bvals
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system

%% Load and visualize the node data
run('points_rhombus_chain.m')

% Determine the adjacency matrix & Total Number of Arches
[data] = determine_adjacency_matrix(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

%% Load and initialize the shape data
data = rise_rhombus_chain_constant(data, bval);

% Start with elastic deformation
[data] = initialize_chiral(zeros(data.N,1),zeros(data.V,1),data);

% Take a look at the initial condition
plot_system_once(data.A0,data)

%% Prepare for Time Integration
data.beta = 0.0075;
T_end = 5000;
data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);
data.A0hat = determine_Ahat_from_A(data.A0,data);

%% Run Time Integration
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

data_Orig = deepCopyStruct(data);
%% Recover A
A_Orig = determine_A_from_Ahat(Ahat,data);
for beta = betavals
    A = A_Orig;
    data = data_Orig;
    %% Visualize the Results
    % data.frames = 100;
    % data.file_name = "b = 0.15 Testing Init";
    % plot_system_over_time(t,A,data)

    %% Look at end of time integration
    plot_system_once(A(end,:)',data)

    data.A0 = A(end,:)';
    data.A0(data.N*data.N_modes+1:end) = 0;

    %% Prepare for Time Integration

    data.beta = beta;

    %data.impose_rotation_at(24) = 1;
    %data.rotation_omega(24) = 2*pi/T_end;
    %data.rotation_mag(24) = 2.5;

    %data.impose_rotation_at(25) = 1;
    %data.rotation_omega(25) = 2*pi/T_end;
    %data.rotation_mag(25) = 2.5;

    data.impose_rotation_at(2) = 1;
    data.rotation_omega(2) = 0.0013;
    data.rotation_mag(2) = 0;

    data.impose_rotation_at(3) = 1;
    data.rotation_omega(3) = 0.0013;
    data.rotation_mag(3) = 0;

    data.impose_rotation_at(46) = 1;
    data.rotation_omega(46) = 0.0013;
    data.rotation_mag(46) = 0;

    data.impose_rotation_at(47) = 1;
    data.rotation_omega(47) = 0.0013;
    data.rotation_mag(47) = 0;



    data.impose_displacement_at(52) = 0.5;
    data.displacement_omega(52) = 2*pi/T_end;

    data = determine_coefficient_matrix(data);
    data = determine_starting_vals(data);
    data = determine_modes_to_skip(data);
    data.A0hat = determine_Ahat_from_A(data.A0,data);

    %% Run Time Integration
    [t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

    %% Recover A
    A = determine_A_from_Ahat(Ahat,data);

    %% Visualize the Results
    data.frames = 100;
    data.file_name = "b = " + num2str(bval) + " beta = " + num2str(data.beta);
    plot_system_over_time(t,A,data)
end
end