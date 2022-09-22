function [shear_rate] = approx_shear_rate(u_old,Hu,flow_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%shear rate (approx. assumes a parabolic profile in vertical)
% shear_rate_alt = 3*u_old./Hu;

shear_rate = (1 + 1/flow_index)/(1 - flow_index/(2*flow_index+1)) * u_old./Hu;

end

