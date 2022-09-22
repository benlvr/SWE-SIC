function [eta_act,u_act] = analytic_soln(h_l,h_r,u_l,u_r,x0,xc,xe,t,g)
%UNTITLED7 Summary of this function goes here
%   Works for when h_r = 0 too.

%middle states
h_m =  1 / (16*g) * (u_l - u_r + 2*sqrt(g*h_l) + 2*sqrt(g*h_r))^2;
u_m = 0;

%celerity of middle state
c_m = sqrt(g*h_m);

%x ranges for pieces of solution
x_a = x0 + t*(u_l - sqrt(g*h_l));
x_b = x0 + t*(u_l + 2*sqrt(g*h_l) - 3*c_m);
x_c = x0 + t*(u_r - 2*sqrt(g*h_r) + 3*c_m);
x_d = x0 + t*(u_r + sqrt(g*h_r));

%rarefections in h
h_raref_1 = 1/(9*g)*(u_l + 2*sqrt(g*h_l) - (xc-x0)/t).^2 .* (xc > x_a) .* (xc < x_b); %connects between h_l and h_m
h_raref_2 = 1/(9*g)*(u_r - 2*sqrt(g*h_r) - (xc-x0)/t).^2 .* (xc > x_c) .* (xc < x_d); %connects between h_r and h_m

%rarefactions in u
u_raref_1 = 2/3*((xe - x0)/t + u_l/2 + sqrt(g*h_l)) .* (xe > x_a) .* (xe < x_b); %connects between u_l and u_m
u_raref_2 = 2/3*((xe - x0)/t + u_r/2 - sqrt(g*h_r)) .* (xe > x_c) .* (xe < x_d); %connects between u_r and u_m

%piecewise solution
eta_act = h_l * (xc < x_a) + h_raref_1 + h_m .* (xc > x_b).*(xc < x_c) + h_raref_2 + h_r * (xc > x_d);
u_act   = u_l * (xe < x_a) + u_raref_1 + u_m .* (xe > x_b).*(xe < x_c) + u_raref_2 + u_r * (xe > x_d);

end

