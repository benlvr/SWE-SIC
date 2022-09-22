function [eta_act,u_act] = act_solution(xc,t,a,C,g,H0,Lx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

eta_act = 0.5*a*cos(pi*xc/Lx)*cos(pi*sqrt(g*H0)*t/Lx);
u_act = 0.5*a*sqrt(g*H0)*(sin(pi*xc/Lx)*sin(pi*sqrt(g*H0)*t/Lx) + C)./ (0.5*a*cos(pi*xc/Lx)*cos(pi*sqrt(g*H0)*t/Lx) + H0);

end

