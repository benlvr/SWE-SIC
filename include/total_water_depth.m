function [H] = total_water_depth(x,eta,wd_tol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%x and eta both cell centered
h = bathymetry(x);
H = h+eta;
% H = max(wd_tol,H);
H = max(0,H);

% p = -1.0/wd_tol^2 * H.^3 + 2.0/wd_tol * H.^2;
% w = find(H > wd_tol);
% d = find(H <= wd_tol);
% F = zeros(size(H));
% F(d) = p(d);
% F(w) = H(w);
% H = F;

%x edge centered and eta cell centered (trap rule)
% h = 0.5*(bathymetry(x(1:end-1)) + bathymetry(x(2:end)));
% H = max(0,h+eta);
end

