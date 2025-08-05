function [data] = determine_time_to_periodic(data)
% Load in time integration left and right
rounding = 4;
points_time_integration = data.points_time_integration;
points = data.points;

up_adjac = triu(data.adjacency_matrix_time_integration,1);
[leftpointnum_t, rightpointnum_t] = find(up_adjac == 1);

% Load in periodic left and right
up_adjac = triu(data.adjacency_matrix,1);
[leftpointnum_p, rightpointnum_p] = find(up_adjac == 1);
leftpoint_p = round(points(leftpointnum_p,:),rounding);
rightpoint_p = round(points(rightpointnum_p,:),rounding);

leftright = [leftpoint_p rightpoint_p];

% (By arch number)
% map_time_to_periodic = [1 3 -2 ... 1]
% 1 time integration -> 1 periodic
% 2 time integration -> 3 periodic
% 3 time integration -> 2 periodic (flipped)
% 56 time integration -> 1 periodic
N_time_integration = size(leftpointnum_t,1);
map_time_to_periodic = zeros(N_time_integration,1);

for time_idx = 1:N_time_integration
    % Find the left and right supports
    lt_point = round(points_time_integration(leftpointnum_t(time_idx,:),:),rounding);
    rt_point = round(points_time_integration(rightpointnum_t(time_idx,:),:),rounding);

    % Shift lt's and rt's x value to match periodic
    while lt_point(1) >= round(data.L_super_cell,rounding)
        lt_point(1) = round(lt_point(1) - round(data.L_super_cell,rounding),rounding);
    end
    while rt_point(1) >= round(data.L_super_cell,rounding)
        rt_point(1) = round(rt_point(1) - round(data.L_super_cell,rounding),rounding);
    end

    % Find the row in periodic time structure that shares the hinges
    % from the time structure

    matched1 = find(all(leftpoint_p == lt_point,2) & all(rightpoint_p == rt_point,2));
    matched2 = -1*find((all(leftpoint_p == rt_point,2) & all(rightpoint_p == lt_point,2)));

    if ~isempty(matched1)
        matched = matched1;
    else
        matched = matched2;
    end
    if isempty(matched)
        disp("Break here");
    end
    % Check if the left and right hinge match
    map_time_to_periodic(time_idx) = matched;

    matched = abs(matched);

    % Check to see if good
    if time_idx > 9999
        f = plot_grid(data,true);
        hold on
        scatter(leftpoint_p(matched,1),leftpoint_p(matched,2),200,'filled','MarkerFaceColor','r')
        scatter(rightpoint_p(matched,1),rightpoint_p(matched,2),200,'filled','MarkerFaceColor','b')

        plot([leftpoint_p(matched,1) rightpoint_p(matched,1)],[leftpoint_p(matched,2) rightpoint_p(matched,2)], ...
            ":",'LineWidth',5,'Color',[252 186 3]/255)

        f_new = figure(2);
        copyobj(allchild(f), f_new);


        % Plot time integration one
        data2.points = data.points_time_integration;
        data2.adjacency_matrix = data.adjacency_matrix_time_integration;
        data2.N = data.N_time_integration;
        data2.N_modes = data.N_modes;
        data2.e_vector = data.e_vector;

        plot_grid(data2,true)
        hold on
        if ~isempty(matched1)
            scatter(lt_point(1),lt_point(2),200,'filled','MarkerFaceColor','r')
            scatter(rt_point(1),rt_point(2),200,'filled','MarkerFaceColor','b')
        else
            scatter(lt_point(1),lt_point(2),200,'filled','MarkerFaceColor','b')
            scatter(rt_point(1),rt_point(2),200,'filled','MarkerFaceColor','r')
        end
        plot([lt_point(1) rt_point(1)],[lt_point(2) rt_point(2)], ...
            ":",'LineWidth',5,'Color',[252 186 3]/255)

        lt_point2 = round(points_time_integration(leftpointnum_t(time_idx,:),:),rounding);
        rt_point2 = round(points_time_integration(rightpointnum_t(time_idx,:),:),rounding);

        if ~isempty(matched1)
            scatter(lt_point2(1),lt_point2(2),200,'filled','MarkerFaceColor','r')
            scatter(rt_point2(1),rt_point2(2),200,'filled','MarkerFaceColor','b')
        else
            scatter(lt_point2(1),lt_point2(2),200,'filled','MarkerFaceColor','b')
            scatter(rt_point2(1),rt_point2(2),200,'filled','MarkerFaceColor','r')
        end
        figure(2)
        disp(time_idx)
    end
end
data.map_time_to_periodic = map_time_to_periodic;
end