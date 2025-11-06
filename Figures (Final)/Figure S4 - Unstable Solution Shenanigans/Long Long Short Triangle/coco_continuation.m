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

%% BP1 
continue_from_BP('vtriangle1',1,'vtriangle1-1',data)


%% BP2
continue_from_BP('vtriangle1',2,'vtriangle1-2',data)

%% BP3
continue_from_BP('vtriangle1',3,'vtriangle1-3',data)


continue_from_BP('vtriangle1-3',1,'vtriangle1-3-1',data)
continue_from_BP('vtriangle1-3',3,'vtriangle1-3-3',data)

%% BP4
continue_from_BP('vtriangle1',4,'vtriangle1-4',data)


%% BP5
continue_from_BP('vtriangle1',5,'vtriangle1-5',data)


continue_from_BP('vtriangle1-5',1,'vtriangle1-5-1',data)
continue_from_BP('vtriangle1-5',3,'vtriangle1-5-3',data)
%% BP6
continue_from_BP('vtriangle1',6,'vtriangle1-6',data)


%% Plot all the Results COCO
idx1 = 1;
idx2 = 2;

figure(9899); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'vtriangle1'        ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vtriangle1-1'      ,'x',idx1,'x',idx2,'b')
coco_plot_bd(theme1, 'vtriangle1-2'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-3'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-4'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-5'      ,'x',idx1,'x',idx2,'b')
% coco_plot_bd(theme1, 'vtriangle1-6'      ,'x',idx1,'x',idx2,'b')


view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

%% Get the relevant data from coco
[a1stab_1,a1unst_1,a2stab_1,a2unst_1,b_1] = get_data_from_coco('vtriangle1');

[a1stab_1_1,a1unst_1_1,a2stab_1_1,a2unst_1_1,b_1_1] = get_data_from_coco('vtriangle1-1');
[a1stab_1_2,a1unst_1_2,a2stab_1_2,a2unst_1_2,b_1_2] = get_data_from_coco('vtriangle1-2');

[a1stab_1_3  ,a1unst_1_3  ,a2stab_1_3  ,a2unst_1_3  ,b_1_3]   = get_data_from_coco('vtriangle1-3');
[a1stab_1_3_1,a1unst_1_3_1,a2stab_1_3_1,a2unst_1_3_1,b_1_3_1] = get_data_from_coco('vtriangle1-3-1');
[a1stab_1_3_2,a1unst_1_3_2,a2stab_1_3_2,a2unst_1_3_2,b_1_3_2] = get_data_from_coco('vtriangle1-3-3');

[a1stab_1_4  ,a1unst_1_4  ,a2stab_1_4  ,a2unst_1_4  ,b_1_4]   = get_data_from_coco('vtriangle1-4');


[a1stab_1_5  ,a1unst_1_5  ,a2stab_1_5  ,a2unst_1_5  ,b_1_5]   = get_data_from_coco('vtriangle1-5');
[a1stab_1_5_1,a1unst_1_5_1,a2stab_1_5_1,a2unst_1_5_1,b_1_5_1] = get_data_from_coco('vtriangle1-5-1');
[a1stab_1_5_2,a1unst_1_5_2,a2stab_1_5_2,a2unst_1_5_2,b_1_5_2] = get_data_from_coco('vtriangle1-5-3');

[a1stab_1_6  ,a1unst_1_6  ,a2stab_1_6  ,a2unst_1_6  ,b_1_6]   = get_data_from_coco('vtriangle1-6');

%%
xlimss = 3*[-0.05 0.05];
ylimss = 3*[-0.05 0.05];
zlimss = 3*[0 0.05];

patchx = [xlimss(1) xlimss(2) xlimss(2) xlimss(1)];
patchy = [ylimss(1) ylimss(2) ylimss(2) ylimss(1)];
patchz = [zlimss(1) zlimss(1) zlimss(2) zlimss(2)];

patch1_color = [50 92 168]/255;
patch2_color = [196 43 43]/255;

patch_alpha = 0.15;

f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2 1.95];

fill3(patchx,[0 0 0 0],patchz,patch1_color,"FaceAlpha",patch_alpha,"EdgeColor","none")
fill3([0 0 0 0],patchy,patchz,patch2_color,"FaceAlpha",patch_alpha,"EdgeColor","none")

plot3(a1stab_1,a2stab_1,b_1,"-","LineWidth",1,"Color","k"); plot3(a1unst_1,a2unst_1,b_1,":","LineWidth",1,"Color","k")

plot3(a1stab_1_2,a2stab_1_2,b_1_2,"-","LineWidth",1,"Color",patch1_color); plot3(a1unst_1_2,a2unst_1_2,b_1_2,":","LineWidth",1,"Color",patch1_color)
plot3(a1stab_1_1,a2stab_1_1,b_1_1,"-","LineWidth",1,"Color",patch2_color); plot3(a1unst_1_1,a2unst_1_1,b_1_1,":","LineWidth",1,"Color",patch2_color)

plot3(a1stab_1_3,a2stab_1_3,b_1_3,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_3,a2unst_1_3,b_1_3,":","LineWidth",1,"Color","k")
plot3(a1stab_1_3_1,a2stab_1_3_1,b_1_3_1,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_3_1,a2unst_1_3_1,b_1_3_1,":","LineWidth",1,"Color","k")
plot3(a1stab_1_3_2,a2stab_1_3_2,b_1_3_2,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_3_2,a2unst_1_3_2,b_1_3_2,":","LineWidth",1,"Color","k")


plot3(a1stab_1_4,a2stab_1_4,b_1_4,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_4,a2unst_1_4,b_1_4,":","LineWidth",1,"Color","k")

plot3(a1stab_1_5,a2stab_1_5,b_1_5,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_5,a2unst_1_5,b_1_5,":","LineWidth",1,"Color","k")
plot3(a1stab_1_5_1,a2stab_1_5_1,b_1_5_1,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_5_1,a2unst_1_5_1,b_1_5_1,":","LineWidth",1,"Color","k")
plot3(a1stab_1_5_2,a2stab_1_5_2,b_1_5_2,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_5_2,a2unst_1_5_2,b_1_5_2,":","LineWidth",1,"Color","k")


plot3(a1stab_1_6,a2stab_1_6,b_1_6,"-","LineWidth",1,"Color","k"); plot3(a1unst_1_6,a2unst_1_6,b_1_6,":","LineWidth",1,"Color","k")


view(3);
grid()
xlim(xlimss); ylim(ylimss); zlim(zlimss)
% xlabel("$a_1$"); ylabel("$a_2$"); zlabel("$a_3$")


set(gca,"FontSize",10)
exportgraphics(gcf,"Long Long Short.png","Resolution",3*600)