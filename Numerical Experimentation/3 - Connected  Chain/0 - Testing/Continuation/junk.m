%% Find expand
data.vertex_map_p2f = [1 13];
data = determine_per_to_finite(data);

%%
A = zeros(data.N*data.N_modes,1);

change_arch_number = 13;
A(change_arch_number*data.N_modes - (data.N_modes-1),1) = 0.75;

plot_system_once(A,data);