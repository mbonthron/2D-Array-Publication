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
% run('points_rhombus.m')
run('points_triangle.m')

% Determine the adjacency matrix & Total Number of Arches
[data] = determine_adjacency_matrix(data);

% Visualize the point and connection between the nodes
plot_grid(data,true);

%% Load and initialize the shape data
% run('rise_chain_decline')
run('rise_triangle_constant')

%% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(N,1),zeros(data.V,1),data);

% Take a look at the initial condition
plot_system_once(data.A0,data)



%% Start for ode45
% data.impose_displacement_at(1) = 0.5;
% data.displacement_omega(1) = 1.2566e-04;

% data.impose_rotation_at(2) = 1;

data = determine_coefficient_matrix(data);

%%
data.A0(1) = 1e-4;
data.beta = 0.005;
[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 50000],data.A0);

plot_system_once(A(end,:),data)


%%
data.impose_rotation_at(1) = 1;
data.rotation_omega(1) = 1.2566e-04; 

% data.impose_rotation_at(3) = 1;
% data.rotation_omega(3) = 1.2566e-04; 

data.A0 = A(end,:)';
data.A0(data.N*data.N_modes+1:end) = 0;

data = determine_coefficient_matrix(data);

data = determine_starting_vals(data);

[t,A] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 4*50000],data.A0);


%%
data.frames = 200;
data.file_name = "Moving 1 Hinge";
plot_system_over_time(t,A,data)

%%
actual = 1*A(:,1)+2*A(:,2)+3*A(:,3);

figure(1);clf ; hold on
plot(t,actual,'LineWidth',2)
plot(t,data.initial_angle(1)*cos(data.rotation_omega(1)*t),":",'LineWidth',2)