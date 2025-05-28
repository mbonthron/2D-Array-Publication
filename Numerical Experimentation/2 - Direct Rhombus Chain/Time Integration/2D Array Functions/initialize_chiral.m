function [data] = initialize_from_flat(impose_displacement_at,impose_rotation_at,data)
%INITIALIZE_CIRAL Defines displacement_vector and force_matrix 
%   a system with plastic deformation such that we start each element in
%   its CHIRAL shape (higher) energy well
%% Load Data
N = data.N;
N_modes = data.N_modes;

%% Determine A0 - starting in Chiral shape
% Load in the Chiral Shape from the periodic System
load("b = 0.15 pi.mat","A0p");
A0 = zeros(2*N*N_modes,1);

%% Go Unit Cell by Unit Cell
unit_cell_count = (N - 1) / 16;
% This relates the non-periodic -> periodic
% nonperiodic_order = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
periodic_order = [1 2 3 4 5 6 7 8 9 11 12 15 16 -10 -13 14];

for i = 0:unit_cell_count-1
    % For the ith unit cell (From the periodic system)
    for j = 1:16
        % For each arch - if first unit cell and first arch
        global_arch_number = 16*i+j;
        
        fprintf("%.0f \t %.0f\n",[periodic_order(j) global_arch_number])
        start_idx = N_modes*global_arch_number - (N_modes - 1);
        end_idx = N_modes*global_arch_number;

        start_idxp = N_modes*abs(periodic_order(j)) - (N_modes - 1);
        end_idxp = N_modes*abs(periodic_order(j));
        
        if periodic_order(j) < 0
            A0(start_idx:end_idx) = -1*A0p(start_idxp:end_idxp);

        else
            A0(start_idx:end_idx) = A0p(start_idxp:end_idxp);
        end
        
    end
end

% Add in the initial shape for the last arch
global_arch_number = N;
start_idx = N_modes*global_arch_number - (N_modes - 1);
end_idx = N_modes*global_arch_number;

start_idxp = N_modes*periodic_order(1) - (N_modes - 1);
end_idxp = N_modes*periodic_order(1);
A0(start_idx:end_idx) = A0p(start_idxp:end_idxp);

data.A0 = A0;

%% Force Vectors and Locations
data.force_magnitude = zeros(N,1);
data.force_eta = zeros(N,1);
data.force_omega = zeros(N,1);

%% Imposed Displacement Values
data.impose_displacement_at = impose_displacement_at;
data.displacement_omega = zeros(N,1);

%% Imposed Rotation Values
data.impose_rotation_at = impose_rotation_at;
data.rotation_omega = zeros(data.V,1);


end

