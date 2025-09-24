clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_single_arch.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1];

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
parameter_names = {'b' 't'};
initial_parameter_values = [0;t_val];

computational_domain = [-0.01 0.125*pi];
UZpoint = [0.025 0.050 0.075 0.1]*pi;

iterations_max = 5000;
h_min = 0.0005;
h_max = 0.001;

% Initial Guess for Continuation
Ahat0 = zeros(2*(data.N*(data.N_modes)-data.constraint_count),1);

%% Run the Initial Continuation Problem
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahat0,parameter_names,initial_parameter_values);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'singlearch1',[],1,parameter_names,computational_domain)


%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 1)
BP = coco_bd_labs('singlearch1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','singlearch1',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'singlearch1-1',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 2)
BP = coco_bd_labs('singlearch1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','singlearch1',BP(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'singlearch1-2',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 3)
BP = coco_bd_labs('singlearch1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','singlearch1',BP(3));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'singlearch1-3',[],1,parameter_names,computational_domain)

%% Plot all the Results COCO
figure(100); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'singlearch1'        ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'singlearch1-1'      ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'singlearch1-2'      ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'singlearch1-3'      ,'x',1,'x',2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

