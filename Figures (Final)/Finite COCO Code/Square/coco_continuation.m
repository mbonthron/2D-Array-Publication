clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_square.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1 1];

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
prob = coco_add_event(prob,'UZ','b',data.UZpoint);

coco(prob,'square0',[],1,data.parameter_names,data.computational_domain)

%% Since we detected HB points - rewrite those as UZ to continue from
add_UZ_to_HB_points('square0',prob,'square1',data)

%% BP 1
continue_from_BP('square1',1,'square1-1',data)

%% BP 2
continue_from_BP('square1',2,'square1-2',data)

add_UZ_to_HB_points('square1-2',prob,'square1-2',data)
continue_from_UZ('square1-2',1,'square1-2-1',data)
continue_from_UZ('square1-2',2,'square1-2-2',data)

%% BP 3
continue_from_BP('square1',3,'square1-3',data)

%% BP4
continue_from_BP('square1',4,'square1-4',data)

add_UZ_to_HB_points('square1-4',prob,'square1-4',data)
continue_from_UZ('square1-4',1,'square1-4-1',data)
continue_from_UZ('square1-4',2,'square1-4-2',data)


%% HB 1
continue_from_UZ('square1',1,'square1-5',data)
continue_from_BP('square1-5',1,'square1-5-1',data)
continue_from_BP('square1-5',2,'square1-5-2',data)

%% HB 2
continue_from_UZ('square1',2,'square1-6',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'square1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-1'        ,'x',idx1,'x',idx2,'b')

coco_plot_bd(theme1, 'square1-2'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-2-1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-2-2'        ,'x',idx1,'x',idx2,'b')

coco_plot_bd(theme1, 'square1-3'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-4'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-4-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-4-2'      ,'x',idx1,'x',idx2,'b')

coco_plot_bd(theme1, 'square1-5'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-5-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'square1-5-2'      ,'x',idx1,'x',idx2,'b')


coco_plot_bd(theme1, 'square1-6'      ,'x',idx1,'x',idx2,'b')



view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

