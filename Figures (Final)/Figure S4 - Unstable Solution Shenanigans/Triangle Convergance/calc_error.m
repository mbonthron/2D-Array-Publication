function [error] = calc_error(x,y,yf)

num = 0;
den = 0;

for i = 1:3
    diff_abs = abs(y(i,:) - yf(i,:));

    num = num + trapz(x,diff_abs);
    den = den + trapz(x,abs(yf(i,:)));

    
end
error = num / den;
end

