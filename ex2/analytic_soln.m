function [eta_act,u_act] = analytic_soln(h_l,h_r,x0,xc,xe,t,g,n_tol)
%UNTITLED7 Summary of this function goes here
%   Works for when h_r = 0 too.

flipped = 0;
if h_l < h_r
    h_l_tmp = h_l;
    h_l = h_r;
    h_r = h_l_tmp;
    
    flipped = 1;
end

%coefficients for sextic polynomial for c_m
c4 = -9*g*h_r;
c3 = 16*g*h_r*sqrt(g*h_l);
c2 = -8*g^2*h_r*h_l - g^2*h_r^2;
c0 = g^3*h_r^3;
c_m = solve_cubic(c0,c2,c3,c4,sqrt(g*(h_l+h_r)/2),n_tol);

%define different regions of actual solution
x_a = x0 - t*sqrt(g*h_l);
x_b = x0 + t*(2*sqrt(g*h_l) - 3*c_m);
if c_m == 0
    x_c = 0;
else
    x_c = x0 + t*(2*c_m^2*(sqrt(g*h_l) - c_m))/(c_m^2 - g*h_r);
end

%parts of analytical solution
h_2 = 4/(9*g)*(sqrt(g*h_l) - (xc-x0)/(2*t)).^2;
h_3 = c_m^2/g;
u_2 = 2/3*((xe - x0)/t + sqrt(g*h_l));
u_3 = 2*(sqrt(g*h_l) - c_m);

u_2 = u_2 .* (xe > x_a) .* (xe <= x_b);
u_3 = u_3 .* (xe > x_b) .* (xe < x_c);
h_2 = h_2 .* (xc > x_a) .* (xc <= x_b);
h_3 = h_3 .* (xc > x_b) .* (xc < x_c);

%assemble different pieces of soln
eta_act = h_l .* (xc <= x_a) + h_r .* (xc >= x_c) + h_2 + h_3;
u_act = u_2 + u_3; %soln in other regions is zero

if flipped == 1
    eta_act = flipud(eta_act);
    u_act = -flipud(u_act);
end

end

