data.variance = ([8.7 8.8 8.7 9.2 8.9]/8.7)';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.2*pi];
data.UZpoint = [data.rise_ND];

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

coco(prob,'vdogbone1',[],1,data.parameter_names,data.computational_domain)

%% BP1 
continue_from_BP('vdogbone1',1,'vdogbone1-1',data)

%%
load("Frustrated State Negative a2.mat")
start_coco_from_another_point(Ahat,b,'vdogbone1-8',data)

%%
load("Frustrated State Positive a2.mat")
start_coco_from_another_point(Ahat,b,'vdogbone1-9',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'vdogbone1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vdogbone1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vdogbone1-8'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vdogbone1-9'      ,'x',idx1,'x',idx2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])
