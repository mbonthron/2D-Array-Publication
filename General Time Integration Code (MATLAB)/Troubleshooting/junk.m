N = data.N;
N_modes = data.N_modes;

% C1 = data.coeff_matrix(end-1,1:N*N_modes);
C2 = data.coeff_matrix(end,1:N*N_modes);

% val1 = C1*A(:,1:N*N_modes)';
val2 = C2*A(:,1:N*N_modes)';

% CV1 = data.coeff_matrix(13,1:N*N_modes)*A(:,1:N*N_modes)';
% CV2 = data.coeff_matrix(14,1:N*N_modes)*A(:,1:N*N_modes)';
% CV3 = data.coeff_matrix(15,1:N*N_modes)*A(:,1:N*N_modes)';

a0 = data.initial_angle;
mag = data.rotation_mag(3);
ome = data.rotation_omega(3);

goal = mag*a0*cos(ome*t);

figure(1); clf; hold on
% plot(t,goal)
plot(t,val2,":","LineWidth",2)

