function [data] = determine_starting_vals(data)
%DETERMINE_STARTING_VALS Summary of this function goes here
%   Detailed explanation goes here
%% Load values from data
e       = data.e;
A0      = data.A0;
N       = 1;
N_modes = data.N_modes;

eta = data.eta;

%% Find Starting height for imposed displacement
initial_height = e.*sin(pi*eta) + diag(sin(pi*eta*(1:1:N_modes))*reshape(A0(1:N*N_modes),[N_modes N]));
data.initial_height = initial_height(eta ~= 0);

end