function [outputArg1,outputArg2] = run_time_integration(shape,data,string_name)
%% Load either the 'S' or 'C' shape
if shape == 'S'
    load('S State.mat','Ahat','b')
    shape_string = "S";
elseif shape == 'C'
    load('C State.mat','Ahat','b')
    shape_string = "C";
else
    fprintf("Gotta pick a shape man")
    quit
end


A = determine_A_from_Ahat(Ahat',data)';
data.A0     = A;
data.A0_D   = data.A0*data.L/pi;

%% Generic Set Up Stuff
data.b_vector = (b*data.variance)';

data.N_cells = 1;
data.plot_grids = 1;

data.arches_to_displace = [1];

data = initialize_time_integration(data);
data.impose_displacement_at(data.arches_to_displace) = data.eta;    % eta value
data.alpha = data.alpha_D;

data = determine_coefficient_matrix_2(data);
data = dimensionalize_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip_FD(data);

%% Find Ahatprime
A0hatprime_D = determine_Ahatprime_from_A_FD(data.A0_D,data);

%% Run Time Integration
T_end = data.mms_to_run / data.alpha_D;

tic
[t,Ahatprime] = ode45(@(t,A) arbitrary_grid_ODE_FD(t,A,data),linspace(0,T_end,data.interp_count),A0hatprime_D);
toc

%%
A = determine_A_from_Ahatprime_FD(Ahatprime', data,t')';

data.frames = 100;
data.file_name = string_name;
close all
% plot_system_over_time(t,A,data)

%% Recover Height and Force information
M_Q = determine_M_Q(t,A,data);
force_idx = data.N_modes*data.N+data.V+data.constraint_count+1;
Q = M_Q(force_idx,:);

displacement = data.alpha*t;
force_rxn    = Q;

energy       = cumtrapz(displacement,force_rxn');

%%
string_name = string_name + ".mat";
save(string_name,"displacement","Q","energy","t","A")

%%
figure(100); clf; hold on; 
plot(displacement*1000,Q,"LineWidth",3)
grid()
% xlim([0 2*data.initial_height*1000])
xlabel("Displacement - [mm]")
ylabel("Force -[N]")
set(gca,'FontSize',18)

figure(2); clf;  hold on
plot(displacement,energy,"LineWidth",3);
% xlim([0 2*data.initial_height])
xlabel("Displacement - [mm]")
ylabel("Energy -[Nm]")
grid()

end

