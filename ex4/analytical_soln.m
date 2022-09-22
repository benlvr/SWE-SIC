function [H_act,dH_dx, u_act] = analytical_soln(Lx,xe,g,g_q)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

H_act = (4/g)^(1.0/3.0)*(1.0 + 0.5*exp(-16.0*(xe/Lx - 0.5).^2));
dH_dx = (4/g)^(1.0/3.0)*(-16.0/Lx*(xe/Lx - 0.5).*exp(-16.0*(xe/Lx - 0.5).^2));
u_act = g_q ./ H_act;

end

