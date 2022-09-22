function [eta_old] = set_eta_gc(eta_old,g1,g2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% eta_old(1) = g1;
eta_old(end) = g2;

% eta_old(1) = eta_old(2);
% eta_old(end) = eta_old(end-1);

end

