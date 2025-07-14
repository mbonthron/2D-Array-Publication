function [A] = determine_A_from_Ahat(time,Ahat,data)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
N_modes         = data.N_modes;
initial_height  = data.initial_height;

eta             = data.eta;
alpha           = data.alpha;

coeff = sin(eta*(2:N_modes)*pi);

a1      = (1/sin(eta*pi))*(initial_height*-alpha*time' - coeff*Ahat(1:(N_modes-1),:));
a1dot   = (1/sin(eta*pi))*(-alpha - coeff*Ahat(N_modes:end,:));

A = [a1 ; Ahat(1:(N_modes-1),:) ; a1dot ; Ahat(N_modes:end,:)];

end