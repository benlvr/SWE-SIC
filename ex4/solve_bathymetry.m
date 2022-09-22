function [B] = solve_bathymetry(H_act,u_act,dH_dx,g,C_f,dx)
%BATHYMETRY Function describing bathymetry
% Solve dB/dx = (1 - u^2/(gh))*dH/dx - S_f

N = length(H_act);

% gamma = g*C_f*u_act.*abs(u_act) ./ H_act.^(1/3);
gamma = manning_friction(H_act,u_act,C_f,g);
RHS = (1 - u_act.^2 ./ (g * H_act)).*dH_dx + gamma .* u_act ./(g*H_act);

B = zeros(N,1);
B(1) = 0;
for i = 1:N-1
    B(i+1) = B(i) + 0.5*dx*(RHS(i+1) + RHS(i));
end

end

