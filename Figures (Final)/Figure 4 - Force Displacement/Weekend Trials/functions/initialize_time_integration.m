function [data] = initialize_time_integration(data)
N = data.N;

%% Force Vectors and Locations
data.force_magnitude = zeros(N,1);
data.force_eta = zeros(N,1);
data.force_omega = zeros(N,1);

%% Imposed Displacement Values
data.impose_displacement_at = zeros(data.N,1);
data.displacement_omega = zeros(N,1);

%% Imposed Rotation Values
data.impose_rotation_at = zeros(data.V,1);
data.rotation_omega = zeros(data.V,1);
data.rotation_mag = ones(data.V,1);

end

