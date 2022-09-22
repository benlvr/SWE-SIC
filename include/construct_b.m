function [b] = construct_b(V,flux_x_ex,flux_x_im,dt,theta)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

b = V(2:end-1) - dt*theta*flux_x_im - dt*(1-theta)*flux_x_ex;

end

