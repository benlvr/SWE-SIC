function [V] = total_vol(x,eta,dx,wd_tol)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%Integrate auxilary porosity fcn over cell
%Assume it's piecewise constant in each cell

H = total_water_depth(x,eta,wd_tol);
V = dx * H;

%trapezoidal rule

end

