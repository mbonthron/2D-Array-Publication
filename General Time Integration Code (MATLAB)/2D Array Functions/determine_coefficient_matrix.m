function [data] = determine_coefficient_matrix(data)
%DETERMINE_COEFFICIENT_MATRIX Returns the coefficient matrix given the
% adjacency matrix, imposed displacement conditions and number of modes
%   INPUTS
%   ===================================================
%   data
%
%   OUTPUTS
%   ===================================================
%   coeff_matrix = coefficient matrix of the system
%% Load Data
adjacency_matrix        = data.adjacency_matrix;
impose_displacement_at  = data.impose_displacement_at;
impose_rotation_at      = data.impose_rotation_at;
N                       = data.N;
N_modes                 = data.N_modes;

%%
% We consider each arch to have a moment on both sides, although some of
% these moments will turn out to be zero
num_moments = N*2;

% Get the number of applied displacements
num_displacements_linear = sum(impose_displacement_at~=0);

% Get the number of applied rotation
num_displacements_rotational = sum(impose_rotation_at~=0);

% Total size of the system
sz = N_modes*N+num_moments+num_displacements_linear+num_displacements_rotational;

% Initialize start of the coeffiecient matrix and right hand side
coeff_matrix = [eye(N_modes*N, sz); zeros(num_moments+num_displacements_linear+num_displacements_rotational,sz)];

% Find connections and treat left as lower in number and right as greater
% in number
up_adjac = triu(adjacency_matrix,1);
[left, right] = find(up_adjac == 1);

% So now, each row in left, right is an arch with it's left and right at a
% vertices in the corresponding graph

%Going to combine into a single vector for indexing purposes
total_moment_vertices = [left;right];


%Now need to determine which Moments equal 0 and which moment is which
%If a vertices only has one 1 in the row, then it's only connected to one
%arch and therefore the moment at that vertices is 0
zero_moment_idx = find(sum(adjacency_matrix,1) == 1);

%Find indexes in the total moment vertices array where they equal 0
total_moment_zero_idx = find(total_moment_vertices==zero_moment_idx);


edges = min(total_moment_vertices) : max(total_moment_vertices)+1;
[counts, values] = histcounts(total_moment_vertices, edges);
Elements = values(counts >= 1);
% NOTE TO SELF MODIFY THIS TO BE A CELL
sum_indexes = {};
for k = 1 : length(Elements)
    sum_indexes{end+1} = [find(total_moment_vertices == Elements(k))]; %all moments at each vertex should add to 0
end

%Now we know where the moments equal 0 and where the moments add to zero,
%so all the pure moment equations
%Note: the first step can be skipped since they also must equal 0, but this
%may be implemented later

%num_zeros = length(zero_moment_idx);

%Populate rows after the identity matrix with semi-identity matrix of zeros
%coeff_matrix(sub2ind([sz, sz],[1:num_zeros]+3*N,zero_moment_idx+3*N)) = 1;

%Populate rows after zero moments with sum of moments
for i=1:length(sum_indexes)
    indices_to_sum = sum_indexes{i};
    for j = indices_to_sum'
        coeff_matrix(N_modes*N + 0 +i, N_modes*N + j) = 1;
        % if j > N %Right moment
        %     coeff_matrix(num_modes*N + 0 +i, num_modes*N + j) = 1;
        % else
        %     coeff_matrix(num_modes*N + 0 +i, num_modes*N + j) = 1;
        % end
    end
end
1; %Code breakpoint for debugging

% Now populating the rest of the equations of motion and the remaining
% constraint equations is next

% Equations of motion
Left_moment_coefficient = -1 * (1:N_modes);
Right_moment_coefficient = (1:N_modes);
Right_moment_coefficient(2:2:end) = -Right_moment_coefficient(2:2:end);


%Go through every arch and every moment
for i=1:N
    %Investigate the polarity here
    row = N_modes*i -(N_modes-1);
    coeff_matrix(row:row+(N_modes-1),N_modes*N+i) = Left_moment_coefficient;
    coeff_matrix(row:row+(N_modes-1),N_modes*N+i+length(left)) = Right_moment_coefficient;
end
1;

% Constraint equations

Left_constraint_coefficient = (1:N_modes);
Right_constraint_coefficient = Right_moment_coefficient;


%Find which arches are connected by the same vertice
connected_arches = {};
%Probably another way to do this but getting the indices of arches
%connected by a vertices

%multi_indexes = sum_indexes(find(counts > 1));

% for i=1:length(multi_indexes)
%     indices = multi_indexes{i};
%     arch_indices =zeros(length(indices),1);
%     for j=1:length(indices)
%         if indices(j) > N
%             arch_indices(j) = indices(j) - N;
%         else
%             arch_indices(j) = indices(j);
%         end
%     end
%     connected_arches{end+1} =  arch_indices;
% end

1;
count = 1;
%number of constraints is number of arches -1
%Need to get pairs of arches in terms of their indices, not their vertexes
%indices

%First need to make arch_ad_matrix with arch indices
%However, the left and thr right matrix already contain the 2D indices of these values

%Then need to see which "rows" including both left and right, contain the
%same values

%Then need some way of identifying left or right for the next step

%Create matrix with left and right so I can later check if each row
%contains a vertex
left_and_right = [left, right];
%Create arches array with numbers corresponding to each arch index
arches = 1:N;
count = 1;
%Create arch pairs array to compile all pairs of arches that intersect at a
%point, this assumes the entire system is connected, however an unconnected
%system would be irrevelant
arch_pairs = zeros(N-1,2);
%Iterate through each vertex
for i = 1:length(adjacency_matrix)
    %Find all arches at a vertex
    all_arches_at_i = arches(any(left_and_right == i,2));
    %If more than 2 arches are at a vertex, continue
    if length(all_arches_at_i) > 1
        left_arches_at_i  = arches(left == i);
        %right_arches_at_i  = arches(right == i);
        j = 1;
        i_count = 1;
        %While j remains an index and more importantly while the number of
        %pairs is less than the total number of pairs/arches at i since the
        %last pair isn't independent with the others
        while j < length(all_arches_at_i) && i_count < length(all_arches_at_i)
            k = j +1;
            if ismember(all_arches_at_i(j),left_arches_at_i)
                %If the arch is connected on left side,
                %change to negative for identification
                %perhaps latter
                new_j = -all_arches_at_i(j);
            else
                new_j = all_arches_at_i(j);
            end
            while k <= length(all_arches_at_i) && i_count < length(all_arches_at_i)
                %left arches will turn negative
                if ismember(all_arches_at_i(k),left_arches_at_i)
                    new_k = -all_arches_at_i(k);
                else
                    new_k = all_arches_at_i(k);
                end
                k = k +1;
                %Save pair of arches at vertex to arch pairs and repeat
                arch_pairs(count, :) = [new_j, new_k];
                count = count +1;
                i_count = i_count+ 1;
            end
            j = j +1;
        end
    end
end

count = 1;
for i=1:size(arch_pairs,1)
    left_count = 0;
    right_count = 0;
    for j=1:2
        if arch_pairs(i,j) < 0
            left_count = left_count +1;
            %turn the index back positive
            col = N_modes* -arch_pairs(i,j) -(N_modes-1);
            %coeff_matrix(num_modes*N+length(sum_indexes)+count,col:col+num_modes-1) = Left_constraint_coefficient;
            if left_count == 1
                coeff_matrix(N_modes*N+length(sum_indexes)+count,col:col+N_modes-1) = Left_constraint_coefficient;
            else %Need to flip the sign of the slope when same
                coeff_matrix(N_modes*N+length(sum_indexes)+count,col:col+N_modes-1) = -Left_constraint_coefficient;
            end
        else
            right_count = right_count +1;
            col = N_modes* arch_pairs(i,j) -(N_modes-1);
            %coeff_matrix(num_modes*N+length(sum_indexes)+count,col:col+num_modes-1) = Right_constraint_coefficient;
            if right_count == 1
                coeff_matrix(N_modes*N+length(sum_indexes)+count,col:col+N_modes-1) = Right_constraint_coefficient;
            else %Need to flip the sign of the slope
                coeff_matrix(N_modes*N+length(sum_indexes)+count,col:col+N_modes-1) = -Right_constraint_coefficient;
            end
        end
    end
    count = count +1;
end
1;

data.constraint_count = size(arch_pairs,1);


%% If there are Imposed Displacements
arches_with_displacements = find(impose_displacement_at~=0);
for i=1:num_displacements_linear
    arch = arches_with_displacements(i);
    frac = abs(impose_displacement_at(arch));
    a_coeffs = zeros(1,N_modes);
    for j=1:N_modes
        a_coeffs(j) = sin(frac*j*pi);
    end

    first_col_row_idx = N_modes*arch-(N_modes-1);

    %Set last rows
    coeff_matrix(N_modes*N+num_moments+i,first_col_row_idx:first_col_row_idx+N_modes-1) = a_coeffs;

    %Set columns on rightmost side of the coeff matrix
    coeff_matrix(first_col_row_idx:first_col_row_idx+N_modes-1,N_modes*N+num_moments+i) = a_coeffs;
end

%% If there are Imposed Rotations
nodes_with_rotation = find(impose_rotation_at~=0);
if ~isempty(nodes_with_rotation)
rotation_count = 1;
    for i = nodes_with_rotation'
        idx = N_modes*N+num_moments+num_displacements_linear + rotation_count;

        %% Find all arches that intersect with the nodes
        all_arches_at_i = arches(any(left_and_right == i,2));
    
        % Add in the unknown moment based on if left or right for each
        % of the arches at the node
        left_arches_at_i = arches(left == i);
        right_arches_at_i = arches(right == i);
    
        % If a left arch add in the left coefficinets
        for j = left_arches_at_i
            row_idx = N_modes*j-(N_modes-1):N_modes*j; 
            coeff_matrix(row_idx,idx) = Left_moment_coefficient';
        end
    
        for j = right_arches_at_i
            row_idx = N_modes*j-(N_modes-1):N_modes*j; 
            coeff_matrix(row_idx,idx) = Right_moment_coefficient';
        end
    
        %% Add in the constraint equation to the bottom row(s) of the
        % coefficient matrix

        % Take the lowest number arch that is at the hinge
        arch_to_add_constraint = min(all_arches_at_i);
        
        % Find the columns that correspond to the points
        col_idx = N_modes*arch_to_add_constraint-(N_modes-1):N_modes*arch_to_add_constraint;
    
        % Rewrite the bottom row based on slope conditions
        if ismember(arch_to_add_constraint,right_arches_at_i)
            coeff_matrix(idx,col_idx) = -1*Right_constraint_coefficient;
        else
            coeff_matrix(idx,col_idx) = Left_constraint_coefficient;
        end
    
        % Add in the moment to the sum of moment equation
        row_idx = N*N_modes + i;
    
        coeff_matrix(row_idx,idx) = 1;
    
        rotation_count = rotation_count + 1;
    end
end

%%
data.coeff_matrix = coeff_matrix;
data.coeff_matrix_modes = coeff_matrix(1:N_modes*N+num_moments,1:N_modes*N+num_moments);

end