function [eta_old] = set_eta_gc(eta_old,dx,t,g_1,g_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%0th order extrapolation
% eta_old(1) = eta_old(2);
% eta_old(end) = eta_old(end-1);

%reflecting
eta_old(1) = -eta_old(2);
eta_old(end) = -eta_old(end-1);

% eta_old(1) = 2*eta_old(2) - eta_old(3); %1st order extrapolation
% eta_old(end) = 2*g_2 - eta_old(end-1);

end

