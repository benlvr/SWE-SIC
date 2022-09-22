function [Gu] = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Gu = u_star(2:end-1) + dt*(nu*Lap_u - g*(1-theta)*d_eta_n_dx);

end

