function [data] = initialize_plastic_deformation(impose_displacement_at,impose_rotation_at,data)
%INITIALIZE_PLASTIC_DEFORMATION Defines displacement_vector and force_matrix 
%   a system with plastic deformation such that we start each element in
%   its HIGHER energy well
%
%   INPUTS
%   ===================================================
%   N: Number of arches
%   N_modes: Number of modes per arch
%   b_vector: Nx1 vector where each row describes b from \deltaL = b^2/4
%   e_vector: Nx1 vector where each row describes e from plastic deformation
%   t_vector: Nx1 vector where each row describes t from thickness
%
%   OUTPUTS
%   ===================================================
%   A0 = 2*N*N_modesx1 vector representing the system in the initially flat
%       position
%   displacement_vector = Nx1 vector represent where each element has a
%       displacement imposed. If the i value is 0 then no displacement is
%       imposed
%   force_vector = Nx1 vector representing the magnitude of forces on each
%       element.
%   eta_vector = Nx1 vector representing the location of forces on each
%       element
%%
N = data.N;
N_modes = data.N_modes;
b_vector = data.b_vector;
e_vector = data.e_vector;
t_vector = data.t_vector;

%% Determine A0
data.A0 = zeros(2*N*N_modes,1);

a1_less_stable = (1/6)*(-9*e_vector-sign(e_vector).*sqrt(9*e_vector.^2 - 12*t_vector.^2));

data.A0(1:N_modes:N*N_modes) = a1_less_stable;


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
data.rotation_mag = ones(data.V,1);

end

