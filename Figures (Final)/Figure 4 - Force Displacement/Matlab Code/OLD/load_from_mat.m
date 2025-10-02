function [displacement,Q] = load_from_mat(matfilename,max_dist)
%LOAD_FROM_MAT Summary of this function goes here
%   Detailed explanation goes here

load(matfilename)

%%
idx = find(displacement>max_dist,1,"first");

displacement = displacement(1:idx);
Q = Q(1:idx);



end

