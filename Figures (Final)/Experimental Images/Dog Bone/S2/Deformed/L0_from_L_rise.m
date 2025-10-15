% clear
L = 100;            % Distance between hinges[mm]


b_vector = [8.1 9.3 10.7 9.9 9.3];
L0_vector  = [];
for i = 1:length(b_vector)
    L0_vector(i) = integral(@(x) sqrt(1+(pi/L*b_vector(i)*cos(pi*x/L)).^2),0,L);

    fprintf("%.2f mm -> L0 = %.3f mm\n",[b_vector(i) L0_vector(i)])
end

L0_vector