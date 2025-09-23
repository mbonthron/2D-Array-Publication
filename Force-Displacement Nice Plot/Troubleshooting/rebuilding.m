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
run('points_rhombus.m')

% Determine the adjacency matrix & Total Number of Arches
[data] = determine_adjacency_matrix(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

%% Load and initialize the shape data
% run('rise_chain_decline')
run('rise_rhombus_constant')

%% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(N,1),zeros(data.V,1),data);

% Start with plastic deformation
% [data] = initialize_plastic_deformation(zeros(data.N,1),zeros(data.V,1),data);


% Take a look at the initial condition
plot_system_once(data.A0,data)



%% Start for ode45
% data.A0(1) = 1e-5;

data.beta = 0.001;

data.impose_displacement_at(1) = 0.5;
data.displacement_omega(1) = (2*pi)/50000;

data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);



[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 50000],data.A0);

%%
plot_system_once(A(end,:),data);

data.A0 = A(end,:)';
%%
% data.A0 = A(end,:)';
% data.A0(data.N*data.N_modes+1:end) = 0;


data.impose_rotation_at(2) = 1;
data.impose_rotation_at(3) = 1;

data.impose_displacement_at(1) = 0.5;
data.displacement_omega(1) = (2*pi)/50000;

data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);


data.force_magnitude = data.force_magnitude + 1;
force = data.force_magnitude;

force(2) = 0;
data.force_magnitude = force;

data.force_eta = linspace(0,1,data.N);
data.force_omega = linspace(0,pi,data.N);

data.impose_displacement_at(1) = 0.5;
data.displacement_omega(1) = (2*pi)/50000;

data = determine_coefficient_matrix(data);

%%
data.beta = 0.001;
[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 50000],data.A0);

plot_system_once(A(end,:),data)


%%
data.frames = 100;
data.file_name = "Rhombus Pinned 2 3";
plot_system_over_time(t,A,data)