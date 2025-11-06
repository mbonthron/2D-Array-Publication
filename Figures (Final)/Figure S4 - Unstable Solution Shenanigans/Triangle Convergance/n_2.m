clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_triangle.m')

data.N_modes = 2;

data.points = points;

data = determine_adjacency_matrix(data);
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

data.plot_grids = true;
plot_grid(data,true)


t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1];
data.e_vector = 0*[1 1 1];

rng('default')
data.mu = 1;
data.sigma = 0.0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';

data.variance = [1 1.5 1.5]';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.125*pi];
data.UZpoint = [0.1]*pi;

data.iterations_max = 5000;
data.h_min = 0.5*0.0005;
data.h_max = 0.5*0.001;

% Initial Guess for Continuation
Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

%% Run the Initial Continuation Problem
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,data.parameter_names,data.initial_parameter_values);
prob = coco_set(prob,'cont','ItMX', data.iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);
    prob = coco_add_event(prob,'UZ','b',data.UZpoint);

coco(prob,'vtriangle1',[],1,data.parameter_names,data.computational_domain)

continue_from_BP('vtriangle1',1,'vtriangle1-1',data)

%%
A = plot_shape_from_COCO('vtriangle1-1',data);
A = A(:,2);
save("n=1.mat","A","data")