function [H] = solve_elevation(u_act,C_f,S_0,g_h,dx) %(u_act,C_f,S_0,g_h,dx,g_u,t_f)
%BATHYMETRY Function describing bathymetry
% Solve dH/dx = S_0 - C_f*u^2/H^(4/3)

wd_tol = 1.0e-10;

g_u = u_act(1);

N = length(u_act);
H = zeros(N,1);
H(1) = g_h;
for i = 1:N-1
    gamma1 = C_f * u_act(i)^2 / H(i)^(4.0/3.0);
    if H(i) <= wd_tol
        gamma1 = 0.0;
    end
    K1 = H(i) + 0.5*dx*(S_0 - gamma1);
    u_half = 0.5*(u_act(i)+u_act(i+1));
    gamma2 = C_f * u_half^2 / K1^(4.0/3.0);
    if K1 <= wd_tol
        gamma2 = 0.0;
    end
    K2 = S_0 - gamma2;
    H(i+1) = H(i) + dx*K2;
end


% idx_start = ceil(g_u*t_f/dx);
% for i = idx_start:-1:2
%     gamma1 = C_f * u_act(i)^2 / H(i)^(4.0/3.0);
%     if H(i) <= wd_tol
%         gamma1 = 0.0;
%     end
%     H(i-1) = H(i) - dx*(S_0 - gamma1);
% end

% H = real(H);
H = max(H,0);

end

