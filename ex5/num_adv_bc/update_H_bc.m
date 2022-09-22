function [h] = update_H_bc(g_u,C_f,S_0,g_h,dt,wd_tol)
%BATHYMETRY Function describing bathymetry
% Solve dH/dx = S_0 - C_f*u^2/H^(4/3)


% g_gamma1 = C_f*g_u^2/g_h^(-4.0/3.0);
% if g_h <= wd_tol
%     g_gamma1 = 0;
% end
% K1 = g_h + 0.5*dt*g_u*(g_gamma1 - S_0);
% g_gamma2 = C_f*g_u^2/K1^(-4.0/3.0);
% if K1 <= wd_tol
%     g_gamma2 = 0;
% end
% K2 = g_u*(g_gamma2 - S_0);
% g_h = g_h + dt*K2;

g_gamma = g_u^2*C_f/g_h^(-4.0/3.0);
if g_h <= wd_tol
    g_gamma = 0;
end
h = g_h + dt*g_u*(g_gamma - S_0);

end

