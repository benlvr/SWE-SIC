function [B,dB_dx] = bathymetry(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Lx = 200.0;
S_0 = 10.0;

% B = S_0 * x;
B = S_0 * (x/Lx).^2;

% dB_dx = S_0 * ones(size(x));
dB_dx = 2*S_0 * x/Lx^2;
end

