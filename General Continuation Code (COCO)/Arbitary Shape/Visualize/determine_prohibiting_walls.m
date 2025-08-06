function [prohibiting_walls] = determine_prohibiting_walls(data)
%DETERMINE_PROHIBITING_WALLS Summary of this function goes here
%   Detailed explanation goes here

%% Load values from data
adjacency_matrix = data.adjacency_matrix_finite;

thetaHigh = data.thetaHigh;
thetaLow  = data.thetaLow;

V = data.V;

% Threshold for which below these the hinges are not considered rotating
rotation_threshold = 6; % Degrees

% Vector which represents the number of the hinges
hinges = 1:V;

% Find which hinges are considered stationary
stationary_hinges = hinges(abs(thetaHigh - thetaLow) < rotation_threshold);

% Remove stationary hinges from adjacency_matrix
adjacency_matrix(stationary_hinges,:) = 0;
adjacency_matrix(:,stationary_hinges) = 0;

% Remove any nodes that are not connected at all
floating_points = find(sum(adjacency_matrix) == 0);

% Remove floating rows and columns from adjacency matrix
adjacency_matrix(floating_points,:) = [];
adjacency_matrix(:,floating_points) = [];

% Check if the new adjacency_matrix is connected
G = graph(adjacency_matrix);
bins = conncomp(G);
binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});

if isscalar(binnodes)
    % Connected Graphc
    prohibiting_walls = 0;
else
    % Unconnected graph
    prohibiting_walls = 1;
end

%% Scatter The Prohobitive Points on Figure
figure(gcf)
scatter(data.points_finite(stationary_hinges,1),data.points_finite(stationary_hinges,2),500,"x", ...
    "MarkerEdgeColor","r","LineWidth",5)

end

