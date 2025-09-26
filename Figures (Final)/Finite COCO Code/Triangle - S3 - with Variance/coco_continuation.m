clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_triangle.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1];
data.e_vector = 0*[1 1 1];

rng('default')
data.mu = 1;
% data.sigma = 0.0183;
data.sigma = 0.0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';

% data.variance = [1 1.05 0.98]';

data.variance = [1.4 .9 1]';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.2*pi];
data.UZpoint = [0.025 0.050 0.075 0.11338]*pi;

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

%% BP1 
continue_from_BP('vtriangle1',1,'vtriangle1-1',data)


%% BP2
continue_from_BP('vtriangle1',2,'vtriangle1-2',data)

% %% BP3
% continue_from_BP('vtriangle1',3,'vtriangle1-3',data)
% 
% %% BP4
% continue_from_BP('vtriangle1',4,'vtriangle1-4',data)
% 
% %% BP5
% continue_from_BP('vtriangle1',5,'vtriangle1-5',data)
% 
% %% BP6
% continue_from_BP('vtriangle1',6,'vtriangle1-6',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'vtriangle1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vtriangle1-1'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-2'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-3'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-4'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-5'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-6'      ,'x',idx1,'x',idx2,'b')


view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

