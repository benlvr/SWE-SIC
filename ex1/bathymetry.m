function [h] = bathymetry(x)
%BATHYMETRY Function describing bathymetry

Lx = 25;
a = 10;
b = -0.05;
c = 0.2;

h = -(b*(x-a).^2 + c) .* (x > 8) .* (x < 12);

%to test subgrid
A = c/10;
p = Lx/2^9;
w = pi/p;
h = h + A*cos(w*x);

end

