alpha_D = 1/1000;
beta    = 1;

eta_vector = [0.75 0.5 0.25];
clc
for i = 1:length(eta_vector)
    dogbone_s1(beta,alpha_D,eta_vector(i))
    dogbone_s2(beta,alpha_D,eta_vector(i))
end
