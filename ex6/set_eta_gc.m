function [eta_old] = set_eta_gc(eta_old,g_eta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

g_1 = -g_eta;
g_2 = g_eta;

%set ghost cells
eta_old(1)   = -eta_old(2) + 2*g_1; %2*eta_old(2) - eta_old(3);
eta_old(end) = -eta_old(end-1) + 2*g_2; %2*g_2 - eta_old(end-1); 

end

