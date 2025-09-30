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
data.sigma = 0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';

data.variance = [1 1 1]';

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

coco(prob,'triangle0',[],1,data.parameter_names,data.computational_domain)

add_UZ_to_HB_points('triangle0',prob,'triangle1',data)

%% BP1
prob = continue_from_BP('triangle1',1,'triangle1-1',data);
add_UZ_to_HB_points('triangle1-1',prob,'triangle1-1',data)


continue_from_UZ('triangle1-1',5,'triangle1-1-1',data)
continue_from_UZ('triangle1-1',6,'triangle1-1-2',data)

%% BP 2
continue_from_BP('triangle1',2,'triangle1-2',data);

%%
run("rotation.m");

bsym  = 0.0236;
Asym1 = [ -0.0045   -0.0020    0.0002   -0.0087   -0.0000    0.0003   -0.0045    0.0020    0.0002      0         0         0         0         0         0         0        0         0]';

Asym2 = [   cw*Asym1(1:9) ; zeros(9,1)];
Asym3 = [cw*cw*Asym1(1:9) ; zeros(9,1)];

Ahatsym1 = determine_Ahat_from_A(Asym1,data);
Ahatsym2 = determine_Ahat_from_A(Asym2,data);
Ahatsym3 = determine_Ahat_from_A(Asym3,data);

% Asymmetric Solution
basym = 1.05*0.0236;
Aasym1 = [basym 0 0 0 0.5*basym 0 -basym 0 0      0         0         0         0         0         0         0        0         0]';

Aasym2 = [   cw*Aasym1(1:9) ; zeros(9,1)];
Aasym3 = [cw*cw*Aasym1(1:9) ; zeros(9,1)];

Ahatasym1 = determine_Ahat_from_A(Aasym1,data);
Ahatasym2 = determine_Ahat_from_A(Aasym2,data);
Ahatasym3 = determine_Ahat_from_A(Aasym2,data);

%%
start_coco_from_another_point(Ahatsym1,bsym,'triangle1-3-1a',data);
start_coco_from_another_point(Ahatsym2,bsym,'triangle1-3-1b',data);
start_coco_from_another_point(Ahatsym3,bsym,'triangle1-3-1c',data);

start_coco_from_another_point(-1*Ahatsym1,bsym,'triangle1-3-1d',data)
start_coco_from_another_point(-1*Ahatsym2,bsym,'triangle1-3-1e',data)
start_coco_from_another_point(-1*Ahatsym3,bsym,'triangle1-3-1f',data)

start_coco_from_another_point(Ahatasym1,basym,'triangle1-3-2a',data)
start_coco_from_another_point(Ahatasym2,basym,'triangle1-3-2b',data)
start_coco_from_another_point(Ahatasym3,basym,'triangle1-3-2c',data)

start_coco_from_another_point(-1*Ahatasym1,basym,'triangle1-3-2d',data)
start_coco_from_another_point(-1*Ahatasym2,basym,'triangle1-3-2e',data)
start_coco_from_another_point(-1*Ahatasym3,basym,'triangle1-3-2f',data)


%% ========================================================================
%  BRANCH FROM FIRST RUN - "HB" Point (No. 2)
continue_from_UZ('triangle1',2,'triangle1-4',data)


%% Plot all the Results COCO
figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP'}});
coco_plot_bd(theme1, 'triangle1'        ,'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'triangle1-1'      ,'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'triangle1-1-1'    ,'x',1,'x',2,'b')
% coco_plot_bd(theme1, 'triangle1-1-2'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-2'      ,'x',1,'x',2,'b')

% coco_plot_bd(theme1, 'triangle1-3'      ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1b'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1c'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1d'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1e'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1f'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-2a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2b'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2c'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2d'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2e'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2f'    ,'x',1,'x',2,'b')

% coco_plot_bd(theme1, 'triangle1-4'      ,'x',1,'x',2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

