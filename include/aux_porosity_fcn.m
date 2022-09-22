function [p] = aux_porosity_fcn(x,z,wd_tol)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%cell centered x
h = bathymetry(x);

%edge centered x
% h = 0.5*(bathymetry(x(1:end-1)) + bathymetry(x(2:end)));

p = (h+z > wd_tol);

end

