function [p] = quadratic_interp(y,xq)
%UNTITLED Summary of this function goes here
%   xq are query points. Assumes points are located as follows: (0,y0),
%   (1,y1), and (2,y2). Assume y is size length(xq) x 3. xq are the number
%   of grid cells away from the grid pooint of interest (e.g. xq=0 is on a
%   grid point, xq=1 is one grid cell upstream etc).


p = 0.5*(1-xq).*(2-xq).*y(:,1) + xq.*(2-xq).*y(:,2) - 0.5*xq.*(1-xq).*y(:,3);

% (0.5*(1-xq).*(2-xq))'
% (xq.*(2-xq))'
% (0.5*xq.*(1-xq))'
end

