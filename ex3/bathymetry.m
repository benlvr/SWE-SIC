function [h,dh_dx] = bathymetry(x)
%BATHYMETRY Function describing bathymetry
%   Here we just use sin(x)*sin(y)

Lx = 25;
a = 10;
b = -0.05;
c = 0.2;

s = -0.005;

x1 = 8;
x2 = 12;

h = -(b*(x-a).^2 + c) .* (x > x1) .* (x < x2);
% h = s*x.*(x - Lx);

dh_dx = -2*b*(x - a) .* (x > 8) .* (x < 12);


% h = -c/5/Lx*(x-Lx);
% dh_dx = -c/5/Lx*ones(size(x));

%to test subgrid
% A = c/20;           %amplitude
% p = 2^5;            %number of periods
% w = 2*pi/Lx*p;      %frequency
% h = h + A*sin(w*x);
% 
% dh_dx = dh_dx + w*A*cos(w*x);

end

