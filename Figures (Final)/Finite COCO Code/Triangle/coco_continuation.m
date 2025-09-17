clear; clc
%% 
addpath("functions")
addpath("visualize")

%% Initialize 'data' for the specific problem
run('initialize_triangle.m')
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

coco(prob,'triangle0',[],1,parameter_names,computational_domain)
clc


%% Since we detected HB points - rewrite those as UZ to continue from
bd = coco_bd_read('triangle0');
HBlbls = coco_bd_labs('triangle0', 'HB');

bcrits = zeros(1,length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end

HB_chiral_anti_chiral = bcrits(1);

prob = coco_add_event(prob,'UZ','b',bcrits);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triangle1',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 1)
BP = coco_bd_labs('triangle1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);

coco(prob,'triangle1-1',[],1,parameter_names,computational_domain)

% % Since we detected HB points - rewrite those as UZ to continue from
bd = coco_bd_read('triangle1-1');
HBlbls = coco_bd_labs('triangle1-1', 'HB');
bcrits = zeros(1,length(HBlbls));
for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end

prob = coco_add_event(prob,'UZ','b',UZpoint);
prob = coco_add_event(prob,'UZ','b',bcrits(1)*1.0001); % If we impose the 'UZ' point exactly, the curve is unstable

coco(prob,'triangle1-1',[],1,parameter_names,computational_domain)
 
% Branch switch from the HB points
UZ = coco_bd_labs('triangle1-1', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1-1',UZ(5));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triangle1-1-1',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1-1',UZ(6));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triangle1-1-2',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 2)
BP = coco_bd_labs('triangle1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1',BP(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triangle1-2',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - "HB" Point (No. 1) - Chiral and Anti Chiral
%  Shapes
UZ = coco_bd_labs('triangle1', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1',UZ(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

% Add a UZ point just above the HB point
prob = coco_add_event(prob,'UZ','b',1.05*HB_chiral_anti_chiral);

coco(prob,'triangle1-3',[],1,parameter_names,computational_domain)


%%
% plot_shape_from_COCO('triangle1-3',data)
run("rotation.m");

% bsym  = 0.1*pi;
% Asym1 = [-0.2493   -0.0860    0.0258   -0.3125    0.0000   -0.0105   -0.2493    0.0860    0.0258      0         0         0         0         0         0         0        0         0]';

% bsym  = 0.025*pi;
% Asym1 = [-0.0574   -0.0208    0.0053   -0.0776    0.0000   -0.0018   -0.0574    0.0208    0.0053      0         0         0         0         0         0         0        0         0]';

% Symmetric Solution
bsym  = 1.05*HB_chiral_anti_chiral;
Asym1 = [ -0.0045   -0.0020    0.0002   -0.0087   -0.0000    0.0003   -0.0045    0.0020    0.0002      0         0         0         0         0         0         0        0         0]';

Asym2 = [   cw*Asym1(1:9) ; zeros(9,1)];
Asym3 = [cw*cw*Asym1(1:9) ; zeros(9,1)];

Ahatsym1 = determine_Ahat_from_A(Asym1,data);
Ahatsym2 = determine_Ahat_from_A(Asym2,data);
Ahatsym3 = determine_Ahat_from_A(Asym3,data);

% Asymmetric Solution
basym = 1.05*HB_chiral_anti_chiral;
Aasym1 = [basym 0 0 0 0.5*basym 0 -basym 0 0      0         0         0         0         0         0         0        0         0]';

Aasym2 = [   cw*Aasym1(1:9) ; zeros(9,1)];
Aasym3 = [cw*cw*Aasym1(1:9) ; zeros(9,1)];

Ahatasym1 = determine_Ahat_from_A(Aasym1,data);
Ahatasym2 = determine_Ahat_from_A(Aasym2,data);
Ahatasym3 = determine_Ahat_from_A(Aasym2,data);

% First Symmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatsym1,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-1a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatsym1,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-1b',[],1,parameter_names,computational_domain)

% Second Symmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatsym2,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-2a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatsym2,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-2b',[],1,parameter_names,computational_domain)

% Third Symmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatsym3,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-3a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatsym3,parameter_names,[bsym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-3b',[],1,parameter_names,computational_domain)

% First Asymmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatasym1,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-4a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatasym1,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-4b',[],1,parameter_names,computational_domain)


% Second Asymmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatasym2,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-5a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatasym2,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-5b',[],1,parameter_names,computational_domain)


% Third Asymmetric Solution
prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Ahatasym3,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-6a',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,-1*Ahatasym3,parameter_names,[basym ; t_val]);
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);
coco(prob,'triangle1-3-6b',[],1,parameter_names,computational_domain)

%% ========================================================================
%  BRANCH FROM FIRST RUN - "HB" Point (No. 2)
UZ = coco_bd_labs('triangle1', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triangle1',UZ(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triangle1-4',[],1,parameter_names,computational_domain)

%% Plot all the Results COCO
figure(100); clf; hold on
theme1 = struct('special', {{'HB','BP','EP'}});
coco_plot_bd(theme1, 'triangle1'        ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-1'      ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-1-1'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-1-2'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-2'      ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3'      ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-1b'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-2a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-2b'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-3a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-3b'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-4a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-4b'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-5a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-5b'    ,'x',1,'x',2,'b')

coco_plot_bd(theme1, 'triangle1-3-6a'    ,'x',1,'x',2,'b')
coco_plot_bd(theme1, 'triangle1-3-6b'    ,'x',1,'x',2,'b')


coco_plot_bd(theme1, 'triangle1-4'      ,'x',1,'x',2,'b')

view(3); grid();
xlim([-0.1 0.1]);
ylim([-0.05 0.05])
zlim([0 0.1])

