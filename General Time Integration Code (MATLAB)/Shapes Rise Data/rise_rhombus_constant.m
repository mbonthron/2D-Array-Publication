N = data.N;      % Number of Arches
b = 0.15*pi;    % Constant Rise to be used

data.b_vector = b*ones(N,1);
data.e_vector = 0*ones(N,1);
data.t_vector = 0.01*pi*ones(N,1);
