clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_hexagon.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1 1 1 1];
data.e_vector = 0*[1 1 1 1 1 1];

rng('default')
data.mu = 1;
data.sigma = 0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';

data.variance = [1.0431    1.1864    1.1916    1.0756    1.0525    1.0000]';
data.variance = 0.25*(data.variance-1)+data.variance;

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.125*pi];
data.UZpoint = [0.025 0.050 0.075 0.1]*pi;

data.iterations_max = 5000;
data.h_min = 0.0005;
data.h_max = 0.001;

% Initial Guess for Continuation
Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

%% Run the Initial Continuation Problem
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,data.parameter_names,data.initial_parameter_values);
prob = coco_set(prob,'cont','ItMX', data.iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);

coco(prob,'hexagon0',[],1,data.parameter_names,data.computational_domain)

%% BP1 
continue_from_BP('hexagon1',1,'hexagon1-1',data)

%% BP2
continue_from_BP('hexagon1',2,'hexagon1-2',data)

%% BP3
continue_from_BP('hexagon1',3,'hexagon1-3',data)

%% BP4
continue_from_BP('hexagon1',4,'hexagon1-4',data)

%% BP5
continue_from_BP('hexagon1',5,'hexagon1-5',data)

%% BP6
continue_from_BP('hexagon1',6,'hexagon1-6',data)

%% BP7
continue_from_BP('hexagon1',7,'hexagon1-7',data)

%% BP8
continue_from_BP('hexagon1',8,'hexagon1-8',data)

%% BP9
continue_from_BP('hexagon1',9,'hexagon1-9',data)

%% BP10
continue_from_BP('hexagon1',10,'hexagon1-10',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'hexagon1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-2'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-3'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-4'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-5'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-6'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-7'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-8'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-9'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'hexagon1-10'     ,'x',idx1,'x',idx2,'b')


view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

