clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_dogbone.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1 1 1 1 1];

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

coco(prob,'dogbone0',[],1,parameter_names,computational_domain)

%% Since we detected HB points - rewrite those as UZ to continue from
bd = coco_bd_read('dogbone0');
HBlbls = coco_bd_labs('dogbone0', 'HB');

bcrits = zeros(1,length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end

prob = coco_add_event(prob,'UZ','b',bcrits);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'dogbone1',[],1,parameter_names,computational_domain)


%% BP 1
run('dogbone1_1.m')

%% BP 2
run('dogbone1_2.m')

%% BP 3
run('dogbone1_3.m')

%% BP4
run('dogbone1_4.m')

%% BP5
run('dogbone1_5.m')

%% BP6
run('dogbone1_6.m')

%% BP 7
run('dogbone1_7.m')

%% HB 1
run('dogbone1_8.m')

%% HB 2
run('square1_6.m')

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'dogbone1'        ,'x',idx1,'x',idx2,'b')

coco_plot_bd(theme1, 'dogbone1-1'      ,'x',idx1,'x',idx2,'b')

% coco_plot_bd(theme1, 'dogbone1-2'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-1'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-1-1'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-1-2'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-1-3'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-1-4'    ,'x',idx1,'x',idx2,'b')

% coco_plot_bd(theme1, 'dogbone1-2-2'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-2-1'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-2-2'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-2-3'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-2-2-4'    ,'x',idx1,'x',idx2,'b')

% coco_plot_bd(theme1, 'dogbone1-3'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-3-1'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-3-2'    ,'x',idx1,'x',idx2,'b')

% coco_plot_bd(theme1, 'dogbone1-4'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-4-1'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-4-2'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-4-3'    ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-4-4'    ,'x',idx1,'x',idx2,'b')

% coco_plot_bd(theme1, 'dogbone1-5'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-6'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'dogbone1-7'      ,'x',idx1,'x',idx2,'b')

coco_plot_bd(theme1, 'dogbone1-8'      ,'x',idx1,'x',idx2,'b')


view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

