N = sum(triu(adjacency_matrix,1) ==1,'all');
b = 0.15*pi;    % Constant Rise to be used

b_vector = b*ones(N,1);
e_vector = 0*ones(N,1);
t_vector = 0.01*pi*ones(N,1);