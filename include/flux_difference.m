function [flux_diff] = flux_difference(flux)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flux_diff = flux(2:end) - flux(1:end-1);

end

