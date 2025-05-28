function [data] = determine_starting_vals(data)
%DETERMINE_STARTING_VALS Summary of this function goes here
%   Detailed explanation goes here
%% Load values from data
e_vector = data.e_vector;
A0 = data.A0;
N = data.N;
N_modes = data.N_modes;

%% Find Starting height for imposed displacement
eta = data.impose_displacement_at;
initial_height = e_vector.*sin(pi*eta) + diag(sin(pi*eta*(1:1:N_modes))*reshape(A0(1:N*N_modes),[N_modes N]));
data.initial_height = initial_height(eta ~= 0);

%% Find Starting Angles for imposed rotation
C = sum(data.impose_rotation_at);   % Find the number of vertecies where rotation is imposed
coeff = data.coeff_matrix(end-(C-1):end,1:N*N_modes);          % Load the last C rows of the coefficient matrix

data.initial_angle = coeff*A0(1:N*N_modes);
end