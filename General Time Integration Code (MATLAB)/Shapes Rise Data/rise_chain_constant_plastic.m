rise = 0.15*pi;    % Constant Rise to be used

data.b_vector = 0*ones(N,1);
data.e_vector = rise*ones(N,1);
data.t_vector = 0.01*pi*ones(N,1);

[data] = fix_polarity(data);