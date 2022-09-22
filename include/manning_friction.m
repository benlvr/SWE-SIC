function [gamma] = manning_friction(H,u,C_f,g)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

gamma = g*C_f.*abs(u) ./ H.^(1.0/3.0);

end

