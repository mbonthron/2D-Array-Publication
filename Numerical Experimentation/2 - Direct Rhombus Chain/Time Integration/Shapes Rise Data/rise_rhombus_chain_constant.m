function data = rise_rhombus_chain_constant(data, bval)

    b = bval*pi;    % Constant Rise to be used
    
    data.b_vector = b*ones(data.N,1);
    data.e_vector = 0*ones(data.N,1);
    data.t_vector = 0.01*pi*ones(data.N,1);
end