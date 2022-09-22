function [x_new] = solve_cubic(c0,c2,c3,c4,c_guess,tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Assume d = 0 for guess (should work for low SS flow rates)
x_old = c_guess;

%if one of the initial states is zero, this product is zero
if c0*c2*c3*c4 == 0
    x_new = 0;
else
    e = 10*tol;
    while e > tol
        f  = x_old^6 + c4*x_old^4 + c3*x_old^3 + c2*x_old^2 + c0;
        df = 6*x_old^5 + 4*c4*x_old^3 + 3*c3*x_old^2 + 2*c2*x_old;
        x_new = x_old - f./df;

        e = abs(x_new - x_old);
        x_old = x_new;
    end
end

end

