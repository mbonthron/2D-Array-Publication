addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_dogbone.m')
t_val = data.t_ND;
data.t_vector = t_val*[1 1 1 1 1];
data.e_vector =     0*[1 1 1 1 1];

data.variance = [1 1 1 1 1]';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.2*pi];
data.UZpoint = [0.1]*pi;

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

coco(prob,'dogbone0',[],1,data.parameter_names,data.computational_domain)

%% Since we detected HB points - rewrite those as UZ to continue from
add_UZ_to_HB_points('dogbone0',prob,'dogbone1',data)

%% BP1 - C Mode in Middle
continue_from_BP('dogbone1',1,'dogbone1-1',data)

%% BP2 - S Mode in Middle
continue_from_BP('dogbone1',2,'dogbone1-2',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'dogbone1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'dogbone1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'dogbone1-2'      ,'x',idx1,'x',idx2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.2])

