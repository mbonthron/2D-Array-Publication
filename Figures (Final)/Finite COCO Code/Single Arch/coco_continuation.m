clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_single_arch.m')
t_val = 0.01*pi;
data.t_vector = t_val*[1];
data.e_vector = 0*[1];

rng('default')
data.mu = 1;
data.sigma = 0;
data.variance = normrnd(data.mu,data.sigma,[1,data.N])';
data.variance = [1 1]';

%% Define the function as the arbitrary grid ODE
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

% Set up for COCO specific things
data.parameter_names = {'b' 't'};
data.initial_parameter_values = [0;t_val];

data.computational_domain = [-0.01 0.2*pi];
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

coco(prob,'single_arch1',[],1,data.parameter_names,data.computational_domain)

%% BP1 
continue_from_BP('single_arch1',1,'single_arch1-1',data)

%% BP2
continue_from_BP('single_arch1',2,'single_arch1-2',data)

%% BP3
continue_from_BP('single_arch1',3,'single_arch1-3',data)

%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'single_arch1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-2'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'single_arch1-3'      ,'x',idx1,'x',idx2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])


%% Make the bifurcation diagram for the figure
% Run 1
bd1 = coco_bd_read('single_arch1');
a_data = bd1(2:end,19);

b_run1  = cellfun(@(x) x(1),bd1(2:end,8));
a1_run1 = cellfun(@(x) x(1),a_data);

eigs = bd1(2:end,15);
lambda1 = cellfun(@(x) max(real(x)),eigs);


% Run 2
bd2 = coco_bd_read('single_arch1-1');
a_data = bd2(2:end,19);

b_run2  = cellfun(@(x) x(1),bd2(2:end,8));
a1_run2 = cellfun(@(x) x(1),a_data);

eigs = bd2(2:end,15);
lambda2 = cellfun(@(x) max(real(x)),eigs);

%
[a1_run1_stab,a1_run1_unst] = separate(a1_run1,lambda1);
[a1_run2_stab,a1_run2_unst] = separate(a1_run2,lambda2);


deltaL_run1 = (b_run1 / 2).^2;
deltaL_run2 = (b_run2 / 2).^2;

%%
L_dimensional = 100;

f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [5/3 1.05];
plot(deltaL_run1*100/pi,a1_run1_stab*100/pi,"k-","DisplayName","Stable")
plot(deltaL_run1*100/pi,a1_run1_unst*100/pi,"r-","DisplayName","Unstable")

plot(deltaL_run2*100/pi,a1_run2_stab*100/pi,"k-","HandleVisibility","off")
plot(deltaL_run2*100/pi,a1_run2_unst*100/pi,"r-","HandleVisibility","off")

% xticks([])
% yticks([])

% legend()

xlabel("$\Delta L$ - [mm]")
ylabel("$b$ - [mm]")
set(gca,'FontSize',10)
