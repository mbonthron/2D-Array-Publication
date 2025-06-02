%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
addpath('2D Array Functions')
addpath('Shapes Point Data')
addpath('Shapes Rise Data')
addpath('Visualize')

%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system

%% Load and visualize the node data
% run('points_chain.m')
run('points_triangle.m')

% Determine the adjacency matrix & Total Number of Arches
[data] = determine_adjacency_matrix(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

%% Load and initialize the shape data
run('rise_square_array_constant')

% Start with elastic deformation
[data] = initialize_from_flat(zeros(data.N,1),zeros(data.V,1),data);

% Take a look at the initial condition
plot_system_once(data.A0,data)

%% Prepare for Time Integration
data.beta = 0.005;

data.A0(1) = 1e-10;

data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);
data.A0hat = determine_Ahat_from_A(data.A0,data);

%% Run Time Integration
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 2*10000],data.A0hat);

%% Recover A
A = determine_A_from_Ahat(Ahat,data);

%% Look at end of time integration
plot_system_once(A(end,:)',data)

data.A0 = A(end,:)';
data.A0(data.N*data.N_modes+1:end) = 0;

%% Impose Rotation
% data.rotation_mag(2) = 1.25;
data.rotation_mag(3) = 3;

% data.impose_rotation_at(2) = 1;
data.impose_rotation_at(3) = 1;

% data.rotation_omega(2) = 2*pi/5000;
data.rotation_omega(3) = 2*pi/5000;

data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data.A0hat = determine_Ahat_from_A(data.A0,data);

%% Run Time Integration
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 4*10000],data.A0hat);

%% Recover A
A = determine_A_from_Ahat(Ahat,data);

%% Visualize the Results
data.frames = 100;
data.file_name = "Rhombus Box - Elastic - From Flat";
plot_system_over_time(t,A,data)