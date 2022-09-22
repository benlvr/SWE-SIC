function [T] = construct_T(Hu,eta_old,theta,dt,dx,g,PHI)
%UNTITLED9
%   Detailed explanation goes here

%We have some options for T based on whether our formulation has non-zero
%or zero viscosity. Here do method (1).
%(1): (No viscosity) We use ghost cells for eta instead of boundary values
%of u when incorporating the u^{n+1} term at the boundary. The vector of
%all 1's is not in the null space of T. BUT if we incorporate 0th order
%extrapolation, we do get 1 in the null space.
%(2): (Viscosity) We use boundary values of u when incorporating the u^{n+1}
%term at the boundary. Vector of 1's is in null space of T.

%friction
N = length(eta_old) - 2;
T = zeros(N,N);
for i = 1:N
    
    if i < N
        T(i,i+1) = Hu(i+1)*PHI(i+1);
        T(i,i) = T(i,i) - Hu(i+1)*PHI(i+1);
    end
    
    if i > 1
        T(i,i-1) = Hu(i)*PHI(i);
        T(i,i) = T(i,i) - Hu(i)*PHI(i);
    end
    
end

% for i = 1:n-2
    %method (1)
    %diagonal entries   
%    T(i,i) = -(Hu(i) + Hu(i+1));
%    if i < n - 2
%        T(i,i+1) = Hu(i+1);
%    end
%    if i > 1
%        T(i,i-1) = Hu(i);
%    end

    %method (2)
%     if i < n - 2
%         T(i,i+1) = Hu(i+1)/dx;
%         T(i,i)   = T(i,i) - Hu(i+1)/dx;
%     end
%     if i > 1
%         T(i,i-1) = Hu(i)/dx;
%         T(i,i)   = T(i,i) -Hu(i)/dx;
%     end
% end

%incorporate ghost cells for method (1) (assuming 0th order extrapolation)
% T(1,1)     = T(1,1)     + Hu(1)*PHI(1);
% T(end,end) = T(end,end) + Hu(end)*PHI(end);

T = -g*dt^2/dx*theta^2 * T;
T = sparse(T);

end




