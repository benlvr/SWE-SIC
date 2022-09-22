function [H] = total_water_depth_db(x,eta,B)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%x and eta both cell centered
h = B(2:end-1);
H = max(0,h+eta);

%x edge centered and eta cell centered (trap rule)
% h = 0.5*(bathymetry(x(1:end-1)) + bathymetry(x(2:end)));
% H = max(0,h+eta);
end

