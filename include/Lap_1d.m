function [L] = Lap_1d(u,dx)
%UNTITLED4 Summary of this function goes here
%   For MAC grids u will be (N+1)x(M) and v will be (N)x(M+1), assuming no
%   ghost cells. This function performs operations on the physical domain.
%   Lap_2d(u) will be (N+1)xM and Lap_2d(v) will be Nx(M+1)
 
% [n,m] = size(u);
% L = zeros(n-2,m);

L = (u(3:end) - 2*u(2:end-1) + u(1:end-2))/dx^2;
end

