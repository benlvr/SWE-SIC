function [u_new] = set_u_bc(u_new,dx,g_q,H_l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

gx_1 = g_q/H_l;

u_new(2)     = gx_1;         %(flux_l < 0)  * (u_new(3))     + (flux_l >= 0)  * (gx_1);
% u_new(end-1) = u_new(end-2); %(flux_r > 0)  * (u_new(end-2)) + (flux_r <= 0) * (u_new(end-2));


% n_l = -1; n_r = 1;
% flux_l = u_new(2)^2     * n_l;
% flux_r = u_new(end-1)^2 * n_r;
% u_new(2)     = gx_1;
% u_new(end-1) = u_new(end-2) + dx*hx_2;
end

