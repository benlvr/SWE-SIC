function [L_oo_e] = compute_Loo_error(u1,u2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

L_oo_e = max(abs(u1 - u2));

end

