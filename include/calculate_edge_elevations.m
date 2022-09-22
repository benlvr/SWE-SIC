function [Hu] = calculate_edge_elevations(u_old,H_old)%H_old
%UNTITLED Calculates H at the u and v edges on *physical domain*
%   Could use upwinding for eta, and then evaluate the bathymetry along the
%   edges.

%determine total elevation at edges via upwinding
a_x_p = (u_old(2:end-1) > 0) + (u_old(2:end-1) == 0)*0.5;
a_x_m = (u_old(2:end-1) < 0) + (u_old(2:end-1) == 0)*0.5;

%Case when Be is sampled at cell centers
Hu = a_x_m.*H_old(2:end) + a_x_p.*H_old(1:end-1);
% Hu = a_x_m.*(eta_old(2:end)  + B(2:end)) + a_x_p.*(eta_old(1:end-1)) + B(1:end-1);

%Case when Be is sampled at edge centers
% Hu = a_x_m.*(eta_old(2:end) + Be(2:end-1)) + a_x_p.*(eta_old(1:end-1)  + Be(2:end-1));
% Hu = a_x_m.*eta_old(2:end) + a_x_p.*eta_old(1:end-1)  + Be(2:end-1);

% Hu = max(H_old(2:end),H_old(1:end-1));
% Hu = max(max(H_old(2:end),H_old(1:end-1)),0);

%Simple averaging
% Hu = (H_old(1:end-1) + H_old(2:end))/2;


end

