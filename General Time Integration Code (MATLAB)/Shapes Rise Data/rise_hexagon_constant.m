rise = 0.15*pi;    % Constant Rise to be used

data.b_vector = rise*ones(data.N,1);
data.e_vector = 0*ones(data.N,1);
data.t_vector = 0.01*pi*ones(data.N,1);