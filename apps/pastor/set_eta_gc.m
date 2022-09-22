function [eta_old] = set_eta_gc(eta_old)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%set ghost cells
eta_old(1)   = eta_old(2);
eta_old(end) = eta_old(end-1);

end

