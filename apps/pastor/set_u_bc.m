function [u_new] = set_u_bc(u_new,dx,g_q,H_l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

gx_1 = 0.0;

u_new(2)     = gx_1;
u_new(end-1) = gx_1;

% u_new(end-1) = u_new(end-2);

end

