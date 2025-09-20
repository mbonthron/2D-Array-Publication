function [V] = calculate_energy(data)
%CALCULATE_ENERGY Summary of this function goes here
%   INPUTS
%   ===================================================
%   A: Vector (or matrix) describing the state space of the system
%       Each row corresponds to a time value 
%       t = 0 -> [a1 a2 a3 ... ]
%       t = 1 -> [a1 a2 a3 ... ]
%
%   N: Number of arches in the system
%   b_vector: Vector describing the b value of each arch from elastic
%       buckling
%   t_vector: Vector describing the thickness of each arch
%   N_modes: Number of modes describing each arch
%
%   OUTPUTS
%   ===================================================
%   V: 1 x N vector which each index corresponds to the potential energy 
%       of the corresponding arch number.
%%
A = data.A0;
N = data.N;
b_vector = data.b_vector;
t_vector = data.t_vector;
N_modes = data.N_modes;

V = zeros(1,N);

for i = 1:N
    % Grab the modes corresponding to the ith arch
    indices = (N_modes*i-(N_modes-1)):(N_modes*i);
    
    sum1 = sum((1:N_modes).^2.*A(indices).^2,2);
    sum2 = sum((1:N_modes).^4.*A(indices).^2,2);

    %
    deltaL = (b_vector(i)/2)^2;
    tval = t_vector(i);

    V(i) = (1/48)*(3*sum1*(sum1-8*deltaL)+2*(sum2)*tval^2) + ((1/144))*(-12*deltaL + tval^2)^2;
end
end