function [data] = initialize_time_integration(data, i)
%INITIALIZE_ELASTIC_DEFORMATION Defines displacement_vector and force_matrix for a
%   elastically deformed system.
%   i is optional argument
%
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
%% Switch the data values for time integration now
N_periodic = data.N;

N = data.N_time_integration;
N_modes = data.N_modes;
N_cells = data.N_cells;
expand  = data.expand;
b_val = data.b;
t_val = data.t;
%% Determine A0
% Load in A0hatp from COCO
load(data.shape_name + " b = "+ num2str(b_val/pi) +" pi t = "+num2str(t_val/pi) + "pi.mat","A0hatp")

% Convert to the A0 system 
if nargin == 1
    A0p = determine_A_from_Ahat(A0hatp',data)';
else
    A0p = zeros(44*6,1);
    A0p(3*i-1) = 0.25*pi;
end

COCO_plot_system_once(A0p,data)


% Initialize Empty Vector
A0 = zeros(2*N*N_modes,1);

for i= 1:N
    periodic_num = data.map_time_to_periodic(i);
    periodic_idx = abs(periodic_num);
    % Grab the modes from A0p
    rows_to_grab = (periodic_idx-1)*N_modes+1:periodic_idx*N_modes;
    modes_p = A0p(rows_to_grab);

    if periodic_num < 0
        modes_p(1:2:N_modes) = -1*modes_p(1:2:N_modes);
    end
    rows_to_write = (i - 1)*N_modes+1:i*N_modes;
    A0(rows_to_write) = modes_p;
end

% Loop over N_cells
% count = 0;
% for i = 1:N_cells
%     for k = 1:N_periodic
%         % Find the finite index from periodic
%         finite_idx = abs(expand(k));
% 
%         % From finite to the time integration system
%         time_idx = finite_idx + N_periodic*(i-1);   % Might be N_finite (?)
% 
%         % Grab the modes from A0p
%         rows_to_grab = (k-1)*N_modes+1:k*N_modes;
%         modes_p = A0p(rows_to_grab);
% 
%         % Check if expand(k) is negative
%         if expand(k) < 0
%             modes_p(1:2:N_modes) = -1*modes_p(1:2:N_modes);
%         end
% 
%         % Write to the time integration A0
%         rows_to_write = (time_idx - 1)*N_modes+1:time_idx*N_modes;
%         A0(rows_to_write) = modes_p;
%         count = count + 1;
%     end
% end
% 
% k = 1;
% i = N_cells + 1;
% while count < N
%     % Find the finite index from periodic
%     finite_idx = abs(expand(k));
% 
%     % From finite to the time integration system
%     time_idx = finite_idx + N_periodic*(i-1);   % Might be N_finite (?)
% 
%     % Grab the modes from A0p
%     rows_to_grab = (k-1)*N_modes+1:k*N_modes;
%     modes_p = A0p(rows_to_grab);
% 
%     % Check if expand(k) is negative
%     if expand(k) < 0
%         modes_p(1:2:N_modes) = -1*modes_p(1:2:N_modes);
%     end
% 
%     % Write to the time integration A0
%     rows_to_write = (time_idx - 1)*N_modes+1:time_idx*N_modes;
%     A0(rows_to_write) = modes_p;
%     count = count + 1;
%     k = k + 1;
% end
    

% Overwrite the adjacency_matrix
data.adjacency_matrix = data.adjacency_matrix_time_integration;
data.points = data.points_time_integration;
data.N = data.N_time_integration;
data.V = size(data.points_time_integration,1);

% Find coeffient matrix and modes to skip for time integration
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

A0hat = determine_Ahat_from_A(A0,data);

% Write the values
data.A0 = A0;
data.A0hat = A0hat;

% Update b_vector
data.b_vector = b_val * ones(N,1);
data.e_vector = 0*ones(N,1);
data.t_vector = data.t_vector(1)*ones(N,1);

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

