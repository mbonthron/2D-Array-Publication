function [f] = plot_grid(adjacency_matrix,points,add_labels)
%PLOT_GRID Plots the adjacency matrix with corresponding points matrix
%   INPUTS
%   ===================================================
%   adjacency_matrix:   Adjacency matrix describing the corresponding points in the point matrix
%   points(:,1) = x position of the nodes
%   points(:,2) = y position of the nodes
%   add_labels  = true / false value if the nodes should be labeled
%
%   OUTPUTS
%   ===================================================
%   f = figure handle the grid was plotted on

%% Plot Styles
grid_color = 0.5*[1 1 1];
grid_linewidth = 1;
grid_linestyle = '-';

node_big_circle_color = 0.75*[1 1 1];
node_big_circle_size  = 200;

node_little_circle_color = 'k';
node_little_circle_size = 50;

text_background_color = [1 1 1 0.25];    % Adds a semitransparent background to make the text easier to read

node_number_font_size = 18;
node_number_color = 'k';

arch_number_font_size = 18;
arch_number_color = [191 2 59]/255;

moment_left_font_size = 8;
moment_left_color = [2 172 191]/255;

moment_right_font_size = 8;
moment_right_color = [191 81 2]/255;

%% Define Commonly Used Variables
x = points(:,1);
y = points(:,2);

N = length(x);

width = max(x) - min(x);
height = max(y) - min(y);

major_axis = max([width height]);

x_center = 0.5*(max(x)-min(x));
y_center = 0.5*(max(y)-min(y));

buffer = 0.25*major_axis;

%% Create Figure
f = figure(9898); clf; hold on
daspect([1 1 1])

% xlim(x_center+.65*[-major_axis major_axis])
% ylim(y_center+.65*[-major_axis major_axis])

xlim(x_center+.5*[-width width]+buffer*[-1 1])
ylim(y_center+.5*[-height height]+buffer*[-1 1])

%% Create Grid
% Connect the corresponding nodes per the adjacency_matrix
for i = 1:N
    for j = i+1:N
        if i ~= j && adjacency_matrix(i,j) == 1
            plot([x(i) x(j)],[y(i) y(j)],grid_linestyle,"LineWidth",grid_linewidth,"Color",grid_color)
        end
    end
end

% Plot the nodes
scatter(x,y,node_big_circle_size, "MarkerFaceColor",node_big_circle_color,"MarkerEdgeColor","k")
scatter(x,y,node_little_circle_size, "MarkerFaceColor",node_little_circle_color,"MarkerEdgeColor","k")


%% Add Text
if add_labels
    offset = .05*major_axis;
    
    % Add the node number label
    text(x-offset,y-offset,string(1:length(x)), "Color",node_number_color, 'FontWeight', 'bold', ...
        'FontSize',node_number_font_size,'HorizontalAlignment','center', ...
        'BackgroundColor',text_background_color)
    
    % Add Arch number label
    up_adjac = triu(adjacency_matrix,1);
    arch_count = sum(up_adjac,'all');
    arches = 1:arch_count;
    
    [left, right] = find(up_adjac == 1);
    
    text((x(left)+x(right))/2,(y(left)+y(right))/2,string(arches), "Color",arch_number_color, "FontWeight", 'bold', ...
        'FontSize',arch_number_font_size,'HorizontalAlignment','center', ...
        'BackgroundColor',text_background_color)
    
    
    % Number of nodes, label each left and right
    for node_num=1:N
        L_arch_num = find(left == node_num);
        R_arch_num = find(right == node_num);
        node_left_right_label = "";
        for i =1:length(L_arch_num)
            L_arch = L_arch_num(i);
            node_left_right_label = node_left_right_label + "L"+num2str(L_arch)+" ";
        end
        for i=1:length(R_arch_num)
            R_arch = R_arch_num(i);
            node_left_right_label = node_left_right_label + "R"+string(R_arch)+" ";
        end
        text(x(node_num)-offset,y(node_num)+offset,node_left_right_label,"Color",moment_right_color, "FontWeight", 'bold', ...
            'FontSize',moment_right_font_size,'HorizontalAlignment','center', ...
            'BackgroundColor',text_background_color)
    end
    
    % Number of nodes, label moment
    for node_num=1:length(adjacency_matrix)
        L_arch_num = find(left == node_num);
        R_arch_num = find(right == node_num);
        node_moment_label = "";
        for i =1:length(L_arch_num)
            L_arch = L_arch_num(i);
            node_moment_label = node_moment_label + "M"+num2str(L_arch)+" ";
        end
        for i=1:length(R_arch_num)
            R_arch = R_arch_num(i);
            node_moment_label = node_moment_label + "M"+string(R_arch + length(left))+" ";
        end
        text(x(node_num)+offset,y(node_num)-offset,node_moment_label,"Color",moment_left_color, "FontWeight", 'bold', ...
            'FontSize',moment_left_font_size,'HorizontalAlignment','center', ...
            'BackgroundColor',text_background_color)
    end
end

end