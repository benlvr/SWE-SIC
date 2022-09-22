function [H_c] = calculate_crit_height(H,u,g)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

H_c = (abs(H.*u) / sqrt(g)) ^ (2.0/3.0);

end

