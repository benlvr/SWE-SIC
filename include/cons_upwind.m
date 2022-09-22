function [u_star] = cons_upwind(u_old,dx,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c = u_old(3:end-2);

u_star = (c > 0) .* (u_old(3:end-2).^2 - u_old(2:end-3).^2)/dx + (c < 0) .* (u_old(4:end-1).^2 - u_old(3:end-2).^2)/dx;

u_star = u_old(3:end-2) - dt*u_star/2;

u_star = [u_old(1); u_old(2); u_star; u_old(end-1); u_old(end)];
end

