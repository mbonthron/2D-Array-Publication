%% Clear Everything so there are no stragglers
clear; clc; close all

%% Add the Paths to the Required Functions
startup
addpath("functions")
addpath("visualize")
addpath("Experimental Images")
addpath("Initialize Functions\")

%% Create Empty Data Structure to be Populated
data = struct();
data.N_modes = 3;   % Number of modes used to describe the system
data.N_cells = 1;
data.plot_grids = 1;

%% Run Continuation to Get Stable Configurations at each b
run('initialize_triangle.m')

% Get points for 2D picture
data.arch_center_distance_mm = 100;
data.number_of_points = 3;
file_name = "01.JPG";

[picture_points_height, picture_points_position] = picture2points(data, file_name);

%%
aN_pic = zeros(data.N,data.N_modes);
for arch_num = 1:size(aN_pic,1)
    % Create the coefficient matrix at each time
    coeff_matrix = zeros(data.number_of_points,data.N_modes);
    for j = 1:data.number_of_points
        coeff_matrix(j,:) = [sin(pi*picture_points_position(arch_num,j)) sin(2*pi*picture_points_position(arch_num,j)) sin(3*pi*picture_points_position(arch_num,j))];
    end
    
    % Use linear algebra to determine what each mode is
    aN_pic(arch_num,:) = coeff_matrix \ -picture_points_height(arch_num,:).';
end

aN_pic_nondim = (aN_pic* pi/data.arch_center_distance_mm)';

A = aN_pic_nondim(:)