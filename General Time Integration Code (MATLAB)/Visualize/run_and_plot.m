function [] = run_and_plot(countidx,damping,T_end,force_sign)
%% Makes a nice video plot
data = struct();

%% Things that might need to change from time to time
shapeNum = 4;
t = 0.1*pi;
b = 0.09*pi; 
sigma = 0.0001;

force_magnitude = force_sign*0.08;

frame_count = 300;

data.timeStr = 'Hexagon Chain - Weekend Run';

%% Create needed data structure
data.N_modes        = 3; 
data.N_cells        = 7;
data.plot_grids     = 1;
data.plot_COCO      = 0;
data.plot_videos    = 1;
data.continue       = 1;
data.b              = b;
data.t              = t;
data.t_vector       = t; % Cringe, but works

data.sigma          = sigma;


%% Load the results (the string array)
load("COCO/"+data.timeStr+"/Sweep Data "+data.timeStr+ " .mat")

archd_idx       = find(strcmp(results(1,:),{'Arch Displacement'}));


%% Initialize the shape (and time integration)
data = init_shape_structural(shapeNum, data, countidx);

% Find arches with negative area to force positive
data.arches_to_force_positive = find(results{countidx+2,archd_idx} < -.3)';
% Find arches with negative area to force positive
data.arches_to_force_negative = find(results{countidx+2,archd_idx} > .3)';

%% Prepare for time integration
data = initialize_time_integration(data);
data.b_vector = data.b_vector.*normrnd(1,sigma,[1,data.N])';
data.t_vector       = t*ones(data.N,1);     % 08/11/2025 Time integration N? or periodic N?

data.beta = .1;    % Initial beta for first time integration

data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);
data.A0hat = determine_Ahat_from_A(data.A0,data);

%% Run the first time integration (relaxing boundaries)
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 200],data.A0hat);

% Look at the end of time integration
% Use the ending configuration for next starting configuration
A = determine_A_from_Ahat(Ahat, data);
data.plot_labels = 1;
%[~] = plot_system_once(A(end,:)',data);
data.plot_labels = 0;

data.A0hat = Ahat(end,:)';
clear Ahat;
data.A0hat(data.N*data.N_modes+1:end) = 0;

%% Determine which arches to force / displace
data.impose_rotation_at(data.nodes_to_hold) = 1;
data.rotation_omega(data.nodes_to_hold) = 0.0013;
data.rotation_mag(data.nodes_to_hold) = 0;

data.impose_displacement_at(data.arches_to_displace) = 0.5;
data.displacement_omega(data.arches_to_displace) = 2*pi/T_end;

data.impose_displacement_at(data.arches_to_displace) = 0.5;
data.displacement_omega(data.arches_to_displace) = 2*pi/T_end;

data.impose_rotation_at(data.nodes_to_rotate) = 1;
data.rotation_omega(data.nodes_to_rotate) = 2*pi/T_end;
data.rotation_mag(data.nodes_to_rotate) = 2.5;

data.force_eta(data.arches_to_force_positive) = .5;
data.force_omega(data.arches_to_force_positive) = 0;
data.force_magnitude(data.arches_to_force_positive) = force_magnitude;

data.force_eta(data.arches_to_force_negative) = .5;
data.force_omega(data.arches_to_force_negative) = 0;
data.force_magnitude(data.arches_to_force_negative) = -force_magnitude;

%% Remake data with new actuation
data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);

%% Set Conditions for final time integration
data.beta = damping;

%% Run Time Integration
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

% Recover A
A = determine_A_from_Ahat(Ahat,data);

%% Visualize the Results
data.frames = frame_count;
tinterp = linspace(t(1),t(end),data.frames);
Ainterp = interp1(t,A,tinterp);

% Check if Folder Exists for the run
if ~exist("Videos\"+ data.timeStr + "\", 'dir')
    mkdir("Videos\"+data.timeStr + "\");
end

plot_system_once(A(1,:),data)

%% Plot the system lowkey nice over time

data.file_name = data.timeStr + "\"+ data.shape_name + " b = " + num2str(b) + " beta = " + num2str(data.beta) + " NumCells = "+ num2str(data.N_cells) + " t = "+num2str(data.t);
%%
plot_system_over_time_lowkey_nice(tinterp,Ainterp,data)
end