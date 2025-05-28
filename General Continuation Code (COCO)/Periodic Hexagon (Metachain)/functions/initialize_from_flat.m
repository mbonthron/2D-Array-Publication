function [A0,displacement_vector,force_vector,eta_vector] = initialize_from_flat(N,N_modes,b_vector,e_vector,t_vector)
%INITIALIZE_FROM_FLAT Defines displacement_vector and force_matrix for an initially flat system
%   Note: This method only makes sense in the context of elastically deformed elements   
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

%% Determine A0
A0 = zeros(2*N*N_modes,1);

%% Determine Displacement_vector
displacement_vector = zeros(N,1);

%% Determine Force Vectors and Locations
force_vector = zeros(N,1);
eta_vector = zeros(N,1);

end

