function [arch_delta] = determine_arch_delta(Ahigh,Alow,data)
%DETERMINE_THETAS Summary of this function goes here
%   Determines the theta's of the hinges
%% Load from data
N  = data.N;    % Number of Arches
N_modes = data.N_modes;     % Number of Modes

% Set up empty arch_delta to be populated
arch_delta = zeros(N,1);

x = linspace(0,pi,100);

% Iterate over each arch
for i = 1:N
    % Find the mode coefficients associated with this arch
    Aparthigh = Ahigh(N_modes*i-(N_modes - 1):N_modes*i)';
    Apartlow  = Alow(N_modes*i-(N_modes - 1):N_modes*i)';
    
    whigh = zeros(size(x));
    wlow  = zeros(size(x));

    % Find the shape of each arch
    for j = 1:N_modes
        whigh = whigh + Aparthigh(j)*sin(j*x);
        wlow  = wlow  + Apartlow(j)*sin(j*x);
    end

    % Determine the area difference between the curves
    area_diff = trapz(x,abs(whigh - wlow));
    
    % Save the area diff into the arch delta array
    arch_delta(i) = area_diff;

end

end

