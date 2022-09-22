function [H] = analytical_soln(x,x0,tau_y,rho,g)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

H = sqrt(2.0 * tau_y / (rho * g) * (x0 - x));

end

