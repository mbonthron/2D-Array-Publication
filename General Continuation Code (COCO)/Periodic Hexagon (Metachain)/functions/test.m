tval = 0.01;
b1 =  0.1;
b2 = 0.15;

a1stab1 = sqrt(12*b1^2/4 - tval^2)/sqrt(3);
a1stab2 = sqrt(12*b2^2/4 - tval^2)/sqrt(3);

N_modes = 3;

p = [b1 b1 b2 b2;
    tval tval tval tval];

A = zeros(6*(N_modes-1),4);
A(1,1) = a1stab1;
A(1,2) = a1stab2;
A(1,3) = a1stab1;
A(1,4) = a1stab2;

arbitrary_grid_ODE(A,p,coeff_matrix,N_modes)