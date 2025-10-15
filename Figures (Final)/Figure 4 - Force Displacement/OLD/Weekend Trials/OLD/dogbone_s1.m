function [] = dogbone_s1(beta,alpha_D,eta)
%DOGBONE_S1 Summary of this function goes here
%   Detailed explanation goes here
%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 1;
data.plot_grids = 1;

b_val  = .087*pi;
t_val  = .00762*pi;
data.arches_to_displace = [1];

%% Run Continuation to Get Stable Configurations at each b
run('initialize_dogbone.m')
run_name = 'dogbone1-1';

data.t_vector = t_val*ones(1,data.N);
data.e_vector = 0*zeros(1,data.N);
data.variance = ([8.7 8.8 8.7 9.2 8.9]/8.1)';

%% Need to pick a starting configuration from COCO
plot_shape_from_COCO(run_name,data)
UZ_instance = 1;

bd = coco_bd_read(run_name);
UZ = coco_bd_labs(run_name, 'UZ');

% Constraint length
C = data.constraint_count;
bcrits      = zeros(1,length(UZ));
Ahat        = zeros(2*(data.N*data.N_modes-C),length(UZ));
stability   = zeros(1,length(UZ));

count = 1;
for k = 1:length(UZ)
    Ahat(:,count)    = coco_bd_val(bd,UZ(k),'x');
    bcrits(:,count)    = coco_bd_val(bd,UZ(k),'b');

    count = count + 1;
end

close all

%% Recover the missing modes from the system
run('physcial_constants.m')
data.V = length(data.points);
data.b = b_val;

data = initialize_time_integration(data);

A = determine_A_from_Ahat(Ahat',data)';

data.A0     = A(:,UZ_instance);
data.A0_D   = data.A0*data.L/pi;


% Indentor Speed
data.impose_displacement_at(data.arches_to_displace) = eta;    % eta value
data.eta = data.impose_displacement_at(data.arches_to_displace);
data.alpha = alpha_D;

data.beta = beta;

data = determine_coefficient_matrix_2(data);
data = dimensionalize_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip_FD(data);

%% Find Ahatprime
A0hatprime_D = determine_Ahatprime_from_A_FD(data.A0_D,data);

A0hatprime_D = A0hatprime_D+1e-10*rand(size(A0hatprime_D));

%% Run Time Integration

T_end = 15;

tic
[t,Ahatprime] = ode45(@(t,A) arbitrary_grid_ODE_FD(t,A,data),linspace(0,T_end,500),A0hatprime_D);
toc

%%
A = determine_A_from_Ahatprime_FD(Ahatprime', data,t')';

%% Recover Height and Force information
M_Q = determine_M_Q(t,A,data);
force_idx = data.N_modes*data.N+data.V+data.constraint_count+1;
Q = M_Q(force_idx,:);

displacement = data.alpha*t;
force_rxn    = Q;

energy       = cumtrapz(displacement,force_rxn');

%%
string_name = "dogbone s1 - beta = "+sprintf("%.2f",beta) + " eta = "+sprintf("%.2f",eta) + " alpha = "+sprintf("%.1f",alpha_D*1000)+ "mms";
save(string_name+".mat","displacement","Q","energy","t","A")


data.frames = 100;
data.file_name = string_name;
close all
plot_system_over_time(t,A,data)

%% Save the data
figure(100); clf; hold on; 
plot(displacement,Q,"LineWidth",3)
grid()
xlim([0 2*data.initial_height])
xlabel("Displacement - [mm]")
ylabel("Force -[N]")
set(gca,'FontSize',18)

figure(2); clf;  hold on
plot(displacement,energy,"LineWidth",3);
xlim([0 2*data.initial_height])
xlabel("Displacement - [mm]")
ylabel("Energy -[Nm]")
grid()
end

