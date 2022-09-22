function [u_old] = set_u_gc(u_old)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

u_old(end) = u_old(end-1); %0th order extrapolation
u_old(1)   = u_old(2);     %0th order extrapolation


end

