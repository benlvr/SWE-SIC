function [u_new] = set_u_bc(u_new,g_u)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

u_new(2) = g_u;

% u_new(end-1) = g_u;

end

