function [res] = Casulli_residual(xc,eta_k,T,b,dx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%evaluate V
V = total_vol(xc,eta_k,dx);

%residual
res = V + T*eta_k - b;

end

