function [data] = init_shape(shapeNum, data)

if shapeNum == 1
    %Rhombus
    %% === Load a points data for the system
    run('points_rhombus_direct')
    
    % Determine the adjacency matrix and number of arches
    data = determine_adjacency_matrix(data);
    data.b_vector = zeros(data.N,1);
    %% Start with elastic deformation
    [data] = initialize_elastic_deformation(zeros(data.N,1),zeros(data.V,1),data);

    %Consider what's actually necessary since this is going into COCO
    data.e_vector = 0*ones(data.N,1);
    data.t_vector = 0.01*pi*ones(data.N,1);
    data.shape_name = 'Rhombus';

elseif shapeNum == 2

elseif shapeNum == 3

elseif shapeNum == 4

elseif shapeNum == 5

elseif shapeNum == 6


end

end

