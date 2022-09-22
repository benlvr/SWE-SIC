function [u_star] = construct_u_star(u_old,u_old_old,dx,dt)
%UNTITLED2 Summary of this function goes here
%   Just glorified upwinding here

N = length(u_old) - 3;

u_star = u_old; 

start_idx = 3   - (u_old(2)   < 0); %if flow is neg at left bdry, decrease idx
end_idx   = N+1 + (u_old(N+2) > 0); %if flow is pos at right bdry, increase idx

% c_p = u_old(3:end-2);
% c_m = c_p;

c_p = u_old(start_idx:end_idx);
c_m = c_p;

%Linear terms
% NL_x = max(c_p,0).*(u_old(3:end-2) - u_old(2:end-3))/dx + min(c_m,0).*(u_old(4:end-1) - u_old(3:end-2))/dx;
% u_star(3:end-2) = u_old(3:end-2) - dt * NL_x;

%if flow is pos, then can do up to u_star(end-1), probably fine to
%extrapolate for ghost cells (which are only needed for diffusion, which is
%most likely for numerical stability)

NL_x = max(c_p,0).*(u_old(start_idx:end_idx) - u_old(start_idx-1:end_idx-1))/dx + min(c_m,0).*(u_old(start_idx+1:end_idx+1) - u_old(start_idx:end_idx))/dx;
u_star(start_idx:end_idx) = u_old(start_idx:end_idx) - dt * NL_x;


end


