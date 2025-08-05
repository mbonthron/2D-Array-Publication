%%
clear; close all


tvals = [.2 8];

for t_idx = 1:length(tvals)
    %% Specify Constraints / Desired Dimensions
    b = 3.5;      % [mm] Dimensional Arch Rise
    L = 75;     % [mm] Dimensional Arch Length
    t = tvals(t_idx);   % [mm] Dimensional Distance from centerline of axis to start of tab
    E = 2e9; % Young's Modulus, Pa
    h = 14; % height of arch

    %% Old Code
    mode = 1;

    %  Calculating Arch Dimensions
    %  SEE ATTACHED PDF FOR DIAGRAM OF RELEVANT DIMENSIONS
    syms x

    %% Calculate ``ideal'' intiial length (if directly to centerline of axis)
    L0 = vpaintegral(sqrt(1+(pi/L*b*cos(pi*x/L))^2),x,0,L);

    %% Conversions
    b=b*10^-3;
    L=L*10^-3;
    t=t*10^-3;
    h=h*10^-3;
    L0 = L0*10^-3;

    %% Calculate q
    switch(mode)
        case 1
            w = b*sin(pi*x/L);
        case 2
            w = 1/2*sqrt(b^2-12)*sin(2*pi*x/L);
        otherwise
            error("Specify Valid Mode");
    end

    A = t*h;
    I = 1/12*h*t^3;

    p = A*E/L*(L0-L-1/2*vpaintegral(diff(w)^2,x,0,L, 'RelTol', 1e-32, 'AbsTol', 0));

    q = E*I*diff(w,4)+p*diff(w,2);

    Rigidity(t_idx) = eval(vpaintegral(E*I*diff(w,4),x,0,L, 'RelTol', 1e-32, 'AbsTol', 0));
    Axial(t_idx) = eval(vpaintegral(p*diff(w,2),x,0,L, 'RelTol', 1e-32, 'AbsTol', 0));

    syms x_prime;
    q_tot(t_idx) = eval(vpaintegral(q,x,0,L, 'RelTol', 1e-32, 'AbsTol', 0));
    V = q_tot(t_idx)/2-int(q,x,0,x_prime);

    V_max = subs(V,x_prime,0);

    tau = 3*V_max/(2*A);
    n(t_idx) = eval(40*10^6/(sqrt(3)*tau));

end

tiledlayout(2,1);
nexttile
hold on;
title("log of Static design factor");
ylabel("log Design Factor");
xlabel("Thickness (mm)");
plot(tvals, log10(n),'LineWidth',3);
hold off;

nexttile

hold on;
title("Bending Rigidity and Shear");
ylabel("Total Force N");
xlabel("Thickness (mm)");
plot(tvals, q_tot, '-','DisplayName','Total q','LineWidth',3);
plot(tvals, Rigidity, '--','DisplayName','Bending Rigidity','LineWidth',3);
plot(tvals, Axial, '--','DisplayName','Axial Force','LineWidth',3);
legend('Location','northwest');
hold off;

nowTime = datetime('now');
timeStr = string(datestr(nowTime, 'yyyy-mm-dd_HH-MM-SS'));

if ~exist("Images/"+timeStr, 'dir')
    mkdir("Images/"+timeStr);
end
print(gcf, "Images/"+timeStr+"/Static Analysis.png", '-dpng', '-r600');

%% New Code
restoredefaultpath
startup
addpath('..\..\General Time Integration Code (MATLAB)\Visualize')
addpath('..\..\General Time Integration Code (MATLAB)\2D Array Functions')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\functions\')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\Visualize\')
addpath('..\..\General Continuation Code (COCO)\Arbitary Shape\')
% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 5;
data.plot_grids = 1;
data.points = [0 0; 0 1];
data = determine_adjacency_matrix(data);
data.b_vector = b*10000;
data.V = 2;
data.N = 1;

% Nodes to hold stationary
data.nodes_to_hold = [];
data.arches_to_displace = [];
data.nodes_to_rotate = [1];
data.arches_to_force_positive = [];
data.arches_to_force_negative = [];

% Start with elastic deformation
[data] = initialize_elastic_deformation(zeros(data.N,1),[1 0],data);

%Consider what's actually necessary since this is going into COCO
data.e_vector = 0*ones(data.N,1);

% Determine the coefficient matrix and number of constraints of the system
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);
data.t_vector = t*ones(data.N,1)*100;

data.beta = .1;
T_end = 200;
data = determine_coefficient_matrix(data);
data = determine_starting_vals(data);
data = determine_modes_to_skip(data);
data.A0 = [b;1e-6;0;0;0;0];
data.A0hat = determine_Ahat_from_A(data.A0,data);
[t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

% Look at end of time integration
A = determine_A_from_Ahat(Ahat, data);