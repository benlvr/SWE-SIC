function [eta_old] = set_eta_gc(eta_old,B,g_h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% g_1 = g_q/u_L;
g_2 = g_h - B(end);

%set ghost cells
eta_old(1)   = eta_old(2); %g_1; %(q_l >= 0) * (eta_old(2))     + (q_l < 0) * (h_o_l);
eta_old(end) = g_2;        %(q_r > 0) * (h_o_r)  + (q_r <= 0) * (g_2);

% eta_old(end) = 2*g_2 - eta_old(end-1); %eta_old(end-1) + dx*h_2;
% eta_old(1) = eta_old(2); %- eta_old(1); %1st order extrapolation
end

