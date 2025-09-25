alpha_D = 1/1000;
beta    = 1;

eta_vector = 1-[0.25 0.30 0.40 0.48 0.50];
clc
for i = 1:length(eta_vector)
    dogbone_s1(beta,alpha_D,eta_vector(i))
    dogbone_s2(beta,alpha_D,eta_vector(i))
end
