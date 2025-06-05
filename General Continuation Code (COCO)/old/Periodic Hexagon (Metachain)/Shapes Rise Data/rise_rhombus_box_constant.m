adjacency_matrix = determine_adjacency_matrix(points);
adjacency_matrix = add_periodicity(adjacency_matrix);

N = sum(triu(adjacency_matrix,1) ==1,'all');

rise = 0.35;    % Constant Rise to be used

b_vector = rise*ones(N,1);
e_vector = 0*ones(N,1);
t_vector = 0.1*pi*ones(N,1);


