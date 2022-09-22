function [P] = free_surf_area(x,eta,dx,wd_tol)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%With x,y and dx,dy we can get quadrature points

%Integrate auxilary porosity fcn over cell
%Assume it's piecewise constant in each cell

P = dx * aux_porosity_fcn(x,eta,wd_tol);

end

