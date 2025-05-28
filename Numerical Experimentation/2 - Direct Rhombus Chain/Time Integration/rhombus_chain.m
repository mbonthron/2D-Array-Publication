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
run('points_rhombus_chain.m')

% Determine the adjacency matrix & Total Number of Arches
[data] = determine_adjacency_matrix(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

%% Load and initialize the shape data
run('rise_rhombus_chain_constant')

% Start with elastic deformation
[data] = initialize_chiral(zeros(data.N,1),zeros(data.V,1),data);

% Take a look at the initial condition
plot_system_once(data.A0,data)


%% Prepare for Time Integration
data.beta = 0.005;

data.impose_rotation_at(2) = 1;
data.rotation_omega(2) = 2*1.2566e-04; 

data.impose_rotation_at(3) = 1;
data.rotation_omega(3) = 2*1.2566e-04; 

data.impose_rotation_at(16) = 1;
data.rotation_omega(16) = 2*1.2566e-04; 

data.impose_rotation_at(17) = 1;
data.rotation_omega(17) = 0; 


data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);

%% Run Time Integration
[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 2*50000],data.A0);

%% Visualize the Results
data.frames = 200;
data.file_name = "b = 0.15";
plot_system_over_time(t,A,data)