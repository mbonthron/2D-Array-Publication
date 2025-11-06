function [aNstable,aNunstable] = separate(aN,stability)
%SEPARATE Summary of this function goes here
%   Detailed explanation goes here
stablepoints = stability < 0;
unstablepoints = stability > 0;

aNstable = aN;
aNunstable = aN;

aNstable(unstablepoints) = nan;
aNunstable(stablepoints) = nan;


end

