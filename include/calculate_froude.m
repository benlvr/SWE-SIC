function [Fr] = calculate_froude(H,u,g)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Fr = abs(u) ./ sqrt(g*H);

end

