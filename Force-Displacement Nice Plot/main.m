%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
restoredefaultpath
startup
% addpath('..\General Time Integration Code (MATLAB)\Visualize')
% addpath('..\General Optimization\')
% addpath('..\General Optimization\Shapes Point Data\')
% addpath('..\General Time Integration Code (MATLAB)\2D Array Functions')
% addpath('..\General Continuation Code (COCO)\Arbitary Shape\functions\')
addpath("functions")
addpath("visualize")
%addpath("data")

%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 1;
data.plot_grids = 1;

b_val  = .1*pi;
t_val = .01*pi;
data.arches_to_displace = [1];
beta = .1;



%% Run Continuation to Get Stable Configurations at each b
% Choose which shape
% shapeNum = 8;
% data.static_only = 1;
% data = init_shape(shapeNum, data);
run('initialize_single_arch.m')
run_name = 'single_arch1-1';


data.t_vector = t_val*ones(1,data.N);


%% Need to pick a starting configuration, maybe from COCO?

plot_shape_from_COCO(run_name,data)
select_branch = 2;

bd = coco_bd_read(run_name);
UZ = coco_bd_labs(run_name, 'UZ');

% Constraint length
C = data.constraint_count;
bcrits      = zeros(1,length(UZ));
Ahat        = zeros(2*(data.N*data.N_modes-C),length(UZ));
stability   = zeros(1,length(UZ));

count = 1;
for k = 1:length(UZ)
    %if coco_bd_val(bd,UZ(k),'b') == b_val
        Ahat(:,count)    = coco_bd_val(bd,UZ(k),'x');
        count = count + 1;
    %end
end

%% Recover the missing modes from the system
data.b = b_val;
data.V = length(data.points);

data = initialize_time_integration(data);

A = determine_A_from_Ahat(Ahat',data)';


data.A0 = A(:,select_branch);
data.A0hat = Ahat(:,select_branch)';

clear A

% Now we have the a values for time integration
T_end = min([2000*sqrt(beta/.1)-400,10000]);

data.impose_displacement_at(data.arches_to_displace) = 0.49;
data.displacement_omega(data.arches_to_displace) = 2*pi/T_end;
data.beta = beta;

% T_end = 2*min([2000*sqrt(beta/.1)-400,10000]);


%% Run Time Integration
data = determine_coefficient_matrix_2(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);

[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);
A = determine_A_from_Ahat(Ahat, data);

data.frames = 100;
data.file_name = "test";
plot_system_over_time(t,A,data)

%% Recover Height and Force information
M_Q = determine_M_Q(t,A,data);
force_idx = data.N_modes*data.N+data.V+data.constraint_count+1;
Q = M_Q(force_idx,:);
figure();
plot(t,Q);



