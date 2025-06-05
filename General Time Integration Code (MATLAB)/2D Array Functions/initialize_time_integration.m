function [data] = initialize_time_integration(b_val,data)
%INITIALIZE_ELASTIC_DEFORMATION Defines displacement_vector and force_matrix for a
%   elastically deformed system.
%
%   INPUTS
%   ===================================================
%   impose_displacement_at: eta values for each arch where displacement
%       should be imposed
%   impose_rotation_at: Vertex numbers where displacement should be imposed
%   data
%
%   OUTPUTS
%   ===================================================
%   data
%% Load data
N = data.N;
N_modes = data.N_modes;

%% Determine A0
data.A0 = zeros(2*N*N_modes,1);
data.A0(1:N_modes:N*N_modes) = data.b_vector;

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

