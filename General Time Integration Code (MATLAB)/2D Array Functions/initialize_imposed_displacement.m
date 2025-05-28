function [A0,displacement_vector,force_vector,eta_vector] = initialize_imposed_rotation(impose_rotation_at,N,N_modes,b_vector,e_vector,t_vector)
%INITIALIZE_ELASTIC_DEFORMATION Defines displacement_vector and force_matrix for a
%   elastically deformed system.
%
%   INPUTS
%   ===================================================
%   impose_displacement_at: eta values for each arch where displacement
%       should be imposed
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

%% Determine A0
A0 = zeros(2*N*N_modes,1);
A0(1:N_modes:N*N_modes) = b_vector;


%% Determine Displacement_vector
displacement_vector = zeros(N,1);

add_displacement_to = find(impose_displacement_at ~= 0);
if ~isempty(add_displacement_to)
    for i = add_displacement_to
        displacement_vector(i) = impose_displacement_at(i);
    end
end

%% Determine Force Vectors and Locations
force_vector = zeros(N,1);
eta_vector = zeros(N,1);

end

