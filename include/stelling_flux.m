function [u_star] = stelling_flux(H_old,Hu,u_old,dx,dt,wd_tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%stelling
H_avg = 0.5 * (H_old(1:end-1) + H_old(2:end));
q = Hu.*u_old(2:end-1);
q_avg = 0.5*(q(1:end-1) + q(2:end));
w = H_avg > wd_tol;

%for non-wet
% c_r = q_avg(1:end-1) ./ H_avg(2:end-1);
% c_l = q_avg(2:end)   ./ H_avg(2:end-1);
% 
% u_star = max(c_r,0) .* (u_old(3:end-2) - u_old(2:end-3))/dx + min(c_l,0) .* (u_old(4:end-1)-u_old(3:end-2))/dx;
% u_star = u_old(3:end-2) - dt * u_star;
% u_star = [u_star(1); u_star(1); u_star; u_star(end); u_star(end)];

%for wet dry cases
q_avg_ext = [q_avg(1); q_avg; q_avg(end)];
c_r_ext = q_avg_ext(1:end-1);
c_l_ext = q_avg_ext(2:end);
c_r_ext(w) = c_r_ext(w) ./ H_avg(w);
c_l_ext(w) = c_l_ext(w) ./ H_avg(w);

u_star_ext = max(c_r_ext,0) .* (u_old(2:end-1) - u_old(1:end-2))/dx + min(c_l_ext,0) .* (u_old(3:end)-u_old(2:end-1))/dx;
u_star_ext = u_old(2:end-1) - dt * u_star_ext;
u_star_ext = [u_star_ext(1); u_star_ext; u_star_ext(end)];

u_star = u_star_ext;

%i = 3
%c_r = q_avg(1) / H_avg(2)
%c_l = q_avg(2) / H_avg(2)
%max(c_r,0) * (u(3) - u(2))/dx + min(c_l,0) * (u(4) - u(3))/dx

%i = N+1 (end-2)
%c_r = q_avg(N-1) / H_avg(N) = q_avg(end-1) / H_avg(end-1)
%c_l = q_avg(N)   / H_avg(N) = q_avg(end)   / H_avg(end-1)
%max(c_r,0) * (u(N+1) - u(N))/dx + min(c_l,0) * (u(N+2) - u(N+1))/dx
%=max(c_r,0) * (u(end-2) - u(end-3))/dx + min(c_l,0) * (u(end-1) - u(end-2))/dx

end

