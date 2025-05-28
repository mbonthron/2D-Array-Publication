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
run('rise_rhombus_constant')

%% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(N,1),zeros(data.V,1),data);

%% Start for ode45
data.A0(1) = 1e-20;

data.beta = 0.001;

data = determine_coefficient_matrix(data);

[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 50000],data.A0);

plot_system_once(A(end,:)',data)
%%
start_row = data.N*data.N_modes + data.V + 1;
end_row = data.N*data.N_modes + 2*data.N; 

C = data.coeff_matrix(start_row:end_row,1:data.N*data.N_modes);
constrain = C*(A(:,1:data.N*data.N_modes)');

figure(2)
plot(t,constrain(3:end,:))

