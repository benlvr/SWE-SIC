function [L_2_e] = compute_L2_error(u1,u2,dx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L_2_e  = norm(u1 - u2) * sqrt(dx);

end

