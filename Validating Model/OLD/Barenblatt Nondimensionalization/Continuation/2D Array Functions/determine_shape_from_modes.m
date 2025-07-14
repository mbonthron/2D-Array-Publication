function [x,w] = determine_shape_from_modes(A,e,horiz_length)
%DETERMINE_SHAPE_FROM_MODES Given the shape of modes for one arch returns
%two vectors x and w(x) describing the deformed shape of the arch
%   INPUTS
%   ===================================================
%   A: Row vector corresponding to the magnitude of each mode;
%   e: Magnitude of plastic deformation in the arch
%   horiz_length: Horizontal distance between the end points
%
%   OUTPUTS
%   ===================================================
%   x = array from zero to horiz_length with 100 elements
%   w = array describing the shape of the arch at the point
%%
x = linspace(0,horiz_length,100);
N_modes = length(A);
w = e*sin(x*pi/horiz_length);
for i = 1:N_modes
    w = w + A(i)*sin(i*x*pi/horiz_length);
end
end

