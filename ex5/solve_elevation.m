function [H] = solve_elevation(u_act,C_f,S_0,g_h,dx)
%BATHYMETRY Function describing bathymetry
% Solve dH/dx = S_0 - C_f*u^2/H^(1/3)

N = length(u_act);
H = zeros(N,1);
H(1) = g_h;
for i = 1:N-1
    K1 = H(i) + 0.5*dx*(S_0 - C_f * u_act(i)^2 / H(i)^(4.0/3.0));
    K2 = S_0 - C_f * u_act(i)^2 / K1^(4.0/3.0);
    H(i+1) = H(i) + dx*K2;
end

end

