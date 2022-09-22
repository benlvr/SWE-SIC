function [H] = analytical_soln(x,x0,H0,B,dB_dx,tau_y,rho,g)
%UNTITLED9 Summary of this function goes here
%   At u=0, gHdH/dx = gHd_z/dx - tau_y/rho
%   Reduces to dH/dx = d_z/dx - tau_y/(rho*g*H)

N = length(x);
dx = x(2) - x(1);
wd_tol = 1.0e-5;

H = zeros(N,1);

% S_0 = dB_dx(1);
%going backwards from front position
%first find index of front
x0_idx = ceil(x0/dx);
H(x0_idx) = H0;
for i = x0_idx:-1:2
    
    %f Euler
%     gamma = tau_y/rho/g/H(i);
%     if H(i) < wd_tol
%         gamma = 0;
%     end
%     RHS = dB_dx(i-1) - gamma;
%     H(i-1) = H(i) - dx*RHS;
    
    %MP
    gamma = tau_y/rho/g/H(i);
    if H(i) < wd_tol
        gamma = 0;
    end
    RHS = dB_dx(i) - gamma;
    H_half = H(i) - 0.5*dx*RHS;
    gamma = tau_y/rho/g/H_half;
    if H_half < wd_tol
        gamma = 0;
    end
    RHS = 0.5*(dB_dx(i) + dB_dx(i-1)) - gamma;
    H(i-1) = H(i) - dx*RHS;
    
end

% H = H / 3;

end

