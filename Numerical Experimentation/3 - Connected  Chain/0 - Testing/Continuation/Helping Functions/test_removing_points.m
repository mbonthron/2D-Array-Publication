%% Load Points Data
yy = sqrt(3)/2;
number_of_elements = 1;

points0 = [0 0;
         0.5 yy;
         1 0;
         1.5 yy;
         2 0;
         2.5 yy;
          3 0;
          3.5 yy
          4 0
          4.5 yy];

points = points0;

for i = 2:number_of_elements
    points_to_add = points0(3:end,:) + [4*(i-1)*ones(8,1) zeros(8,1)];
    points = [points ; points_to_add];
end

% run('points_triangle.m')
% run('points_square2.m')
% run('points_rhombus.m')
% run('points_square2.m')
% run('points_unitcell.m')
% run('points_square_chain.m')


data.points = points;
data.V = length(points);


%% Determine Adjacency Matrix
data = determine_adjacency_matrix(data);

%% Remove desired connections
middles = [1 2;
           2 3;
           3 4;
           4 5;
           5 6;
           6 7;
           7 8;
           8 9;
           9 10];

% remove_middles = [1 3 5 7 9];
remove_middles = [3 6 7];
for i = remove_middles
    data = remove_connection(middles(i,1),middles(i,2),data);
end


data.N = sum(triu(data.adjacency_matrix,1) ==1,'all');
data.N_modes = 1;
data.impose_displacement_at = zeros(data.N,1);
data.impose_rotation_at = zeros(data.V,1);

%% Find the coefficient matrix
data = determine_coefficient_matrix(data);

%% Visualize Points
plot_grid(data,true)
axis off

%% Check if the system is frustrated
cc = data.constraint_count;
C = data.coeff_matrix;
last_N_rows = C(end-cc+1:end,1:data.N*data.N_modes);

first_mode_coeff = last_N_rows(:,1:data.N_modes:end);

r = rank(first_mode_coeff);
[m,n] = size(first_mode_coeff);

if r == n
    fprintf("SYSTEM IS FRUSTRATED\n")
else
    fprintf("SYSTEM IS NOT FRUSTRATED\n")
end

%%
x = 0:1:10;
x2 = 0.5:1:10.5;


figure(1); clf; hold on
scatter(x,zeros(size(x)),50,"filled")
scatter(x2,sqrt(3)/2*ones(size(x2)),50,"filled")
axis equal