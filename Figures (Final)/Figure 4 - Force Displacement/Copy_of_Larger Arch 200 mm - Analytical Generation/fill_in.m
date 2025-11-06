function [f1] = fill_in(displacement,force_mean,force_low,force_high,shape)
%FILL_IN Summary of this function goes here
%   Detailed explanation goes here
if isstring(shape) && shape(1) == "C"
    color = [52 93 168]/255;
elseif isstring(shape) && shape(1) == "S"
    color = [196 43 43]/255;
else
    color = shape;
end


f1 = fill([displacement; flipud(displacement)], ...
     [force_high+force_mean; flipud(force_mean-force_low)], ...
     color, ...           
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4,"HandleVisibility","off");
end

