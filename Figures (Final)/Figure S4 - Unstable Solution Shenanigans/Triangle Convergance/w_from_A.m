function [w,x] = w_from_A(A,N)
%%
x = linspace(0,pi,500);
w = zeros(3,500);

for k = 1:3
    for i = 1:N
        disp(N*(k-1)+i)
        w(k,:) = w(k,:) + A(N*(k-1)+i)*sin(i*x);
    end
end
%%
end

