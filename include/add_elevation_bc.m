function [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,bdry)
%ADD_ELEVATION_BC Chnages nonlinear system to introduce an elevation bc
%   Method assumes velocity boundary conditions are supplied on both
%   boundaries. This function changes corresponding entry of T and b so
%   that a boundary condtion of eta(x_bdry) = g_eta is applied


N = length(eta_old)-2;

%determine cell index: bdry=1 => idx=1, bdyr=2 => idx=N
c_idx = (N-1)*(bdry - 1) + 1;

%determine edge index: idx=1 and idx=N+1
e_idx = N*(bdry - 1) + 1;

%change entry in T
T(c_idx,c_idx) = T(c_idx,c_idx) + g*dt^2*theta^2*Hu(e_idx)*PHI(e_idx);

%add entry to the RHS of the nonlinear system
b(c_idx) = b(c_idx) + g*dt^2*theta^2*Hu(e_idx)*PHI(e_idx)*eta_old(c_idx);

end

