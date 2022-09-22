function [u_new] = set_u_bc(u_new,dx,g_q,H_l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if H_l == 0
    c = 0;
else
    c = 1/H_l;
end
gx_1 = g_q * c;
% gx_1 = g_q / H_l;

u_new(2)     = gx_1;
% u_new(end-1) = u_new(end-2);

end

