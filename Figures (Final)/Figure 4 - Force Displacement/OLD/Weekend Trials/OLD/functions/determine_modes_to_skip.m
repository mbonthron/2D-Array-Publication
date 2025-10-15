function [data] = determine_modes_to_skip(data)
%%
coeff_matrix_modes = data.coeff_matrix_modes;
C = data.constraint_count;
N = data.N;
N_modes = data.N_modes;

%%
last_C_rows = coeff_matrix_modes(end-(C-1):end,1:N*N_modes); 

[~,modes_to_skip] = rref(last_C_rows(:,3:end));
data.modes_to_skip = modes_to_skip + 2;

end