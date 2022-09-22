function [eta_k1, r, iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter)
%UNTITLED8 Summary of this function goes here
%   T and b are dependent only on the values of eta at the previous time
%   step, not the current one.

N = length(eta_old(2:end-1));

guess = eta_old;

B = bathymetry(xc);

iter = 0;
eta_k = guess(2:end-1);
eta_k1 = zeros(N,1);

r = 10 * n_tol; %norm(res);
while r > n_tol && iter < max_iter
    
    %evaluate V
    V = total_vol(xc,eta_k,dx,wd_tol);
        
    %evaluate P
    p = free_surf_area(xc,eta_k,dx,wd_tol);
    P = spdiags(p,0,N,N);
    
    %omit rows/colums of zeros
    %a is boolean array for if a cell or at least its neighbors are nonzero
    A = abs(P+T) > wd_tol;
    a = A*ones(N,1); 
    c = find(a > wd_tol);
    d = find(a <= wd_tol);
    
    %
%     di = find(a(2:end-1) <= wd_tol)+ 1;
%     db = find(a(1) <= wd_tol);
%     db = [db find(a(N) <= wd_tol)+N-1];
   
    %update eta (don't update cells that are empty and have empty
    %neighbors)
    res = V + T*eta_k - b;
    d_eta_k1 = (P(c,c)+T(c,c))\(res(c));
    eta_k1(c) = eta_k(c) - d_eta_k1;  
%     eta_k1(d) = eta_k(d);
    eta_k1(d) = b(d)/dx - B(d);

%     eta_k1(di) = eta_k(di);
%     eta_k1(db) = b(db)/dx - B(db);

    
    %residual
    res = V + T*eta_k - b;

    iter = iter + 1;
    r = norm(res);
    eta_k = eta_k1;
end

end
