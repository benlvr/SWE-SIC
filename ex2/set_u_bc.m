function [u_new] = set_u_bc(u_new,dx,u_l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% hx_2 = 0;
gx_1 = u_l;
% gx_2 = 0;

% u_new(2)     = gx_1;
% u_new(end-1) = u_new(end-2) + dx*hx_2;

%outflow
% u_new(2) = u_new(3);
% u_new(end-1) = u_new(end-2);

%wall
% u_new(2) = 0;
% u_new(end-1) = 0;



end

