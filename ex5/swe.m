
%Based on Hunter 2005 paper
%Diffusive wave equation with friction (solves full SWE)

%VARIABLES
%H:   total elevatation of water, H = eta + B
%eta: water elevation above the plane y = 0
%B:   depth of terrain below the plane y = 0
%u:   depth-averaged velocity in x-direction

%EQUATIONS
%Mass: H_t + (H*u)_x = 0
%Momentum: H*(u_t + uu_x) = -g*H*eta_x - gamma*u

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')

clear all
close all

%parameters
t_f = 3600.0;
Lx = 5000.0; %m
g = 9.81; %m/s^2
S_0 = 0.0; %10^-3, bed slope
C_f = (0.01)^2; %m^(-1/3), Manning friction coefficient, C_f = n^2

%boundary and initial data
u0 = 0.0; %m/s
g_u = 1.0;
g_h = 0; %H(0,0), not BC for all times
C = 3.0/7.0 * (g_h)^(7.0/3.0); %const of int for act soln, should be > C_f*u0^3*Lx
bdry = 2;

%Tolerances
n_tol  = 1.0e-8;
wd_tol = 1.0e-6;
max_iter = 10;

%timestepping parameters
u_max = g_u + sqrt(g*g_h);
CFL = 0.125; %stelling up to CFL=0.15, upwind up to CFL=0.5 (small osc's)
theta = 0.5;

%number of cells in each direction
N = 2^8;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = 100;
stride = ceil(t_incr/dt);

nu = 0.05*dx^2; %numerical viscosity

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize total water column, bathymetry, and surface elevation
H = (7.0/3.0*(C - C_f*u0^3*xc)).^(3.0/7.0);
B = bathymetry(xc);
eta = H - B;
      
%axis
eta_axis = [-2*dx Lx+2*dx -0.5 1.5];
u_axis   = [-2*dx Lx+2*dx -0.5*g_u 1.5*g_u];
q_axis   = [-2*dx Lx+2*dx -0.5*g_u 1.5*g_u];

%initialize bathymetry
B = [B(1); B; B(end)];

%initialize variables
u = u0 * ones(N+1,1);

%create ghost cell arrays
eta_new = [eta(1); eta; eta(end)];
H_new = eta_new + B;
u_new = [g_u; u; 0.0];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

u_old_old = u_old;

gamma = zeros(N+1,1);
Hu = calculate_edge_elevations(u_old,H_old);

t_idx = 0;
t = 0;
while t <= t_f
    
    %set ghost cells
    g_h = (7.0/3.0*(C + C_f*g_u^4*t))^(3.0/7.0);
    g_q = g_u*g_h;
    eta_old = set_eta_gc(eta_old);
    u_old = set_u_gc(u_old);
    
    %plot
    plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
    u_star = stelling_flux(H_old,Hu,u_old,dx,dt,wd_tol);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_star,dx);

    %Casulli's "G" operator
    Gu = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_star,H_old);
      
    %wet and dry vectors
    w = find(Hu > wd_tol);  
    
    %friction
    gamma    = zeros(N+1,1);
    gamma(w) = manning_friction(Hu(w),u_old(w+1),C_f,g);
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);
    
    %construct fluxes
    [flux_ex,flux_im] = construct_lr_fluxes(Hu,Hu,u_old,Gu,PHI);
    [flux_ex,flux_im] = set_flux_bc_alt(flux_ex, flux_im,g_q);
    [flux_diff_ex] = flux_difference(flux_ex);
    [flux_diff_im] = flux_difference(flux_im);
        
    %update volume
    V_new = V_old(2:end-1) - dt*flux_diff_ex;
    
    %reconstruct total elevation H
    H_new(2:end-1) = V_new / dx;
    
    %calculate surface elevation eta
    eta_new(2:N+1) = H_new(2:end-1) - B(2:end-1);
    
    %nonlinear
    T = construct_T(Hu,eta_old,theta,dt,dx,g,PHI);
    b = construct_b(V_old,flux_diff_ex,flux_diff_im,dt,theta);
    [T,b] = add_elevation_bc(Hu,eta_old,PHI,T,b,g,dt,theta,bdry);
    [eta_update,r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);
    eta_new(2:end-1) = eta_update;
    
    eta_new = set_eta_gc(eta_new);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
        
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;
    u_new(2:end-1) = PHI.*RHS;
      
    %set u bc's
    u_new = set_u_bc(u_new,g_u);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update time
    t_idx = t_idx + 1;
    t     = t + dt;
        
    %update previous values to new values
    u_old_old = u_old;
    eta_old = eta_new;
    u_old = u_new;
    
end

%analytical solution
u_act_c = g_u*ones(N,1) .* (xc < g_u*t);
H_act = (7.0/3.0*(C - C_f*u_act_c.^3.*(xc - u_act_c*t))).^(3.0/7.0);
eta_act = H_act - B(2:end-1);
u_act = g_u*ones(N+1,1) .* (xe < g_u*t);

%error in SS soln
L_2_e  = [norm(eta_act - eta_new(2:end-1)) norm(u_act - u_new(2:end-1))] * sqrt(dx)
L_oo_e = [max(abs(eta_act - eta_new(2:end-1))) max(abs(u_act - u_new(2:end-1)))]

%plot final state
plot_state(xc,xe,t,0,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)

%plot for export
figure(2)
lw = 4;
fs = 24;
subplot(1,3,1)
hold on
plot(xc,eta_new(2:end-1),'LineWidth',lw)
plot(xc,eta_act,'--','LineWidth',lw)
plot(xc,-B(2:end-1),'--','LineWidth',lw)
xlabel('x direction')
ylabel('\eta')
legend('numerical','actual','bathymetry','FontSize',fs)
subplot(1,3,2)
hold on
plot(xe, u_new(2:end-1),'LineWidth',lw)
plot(xe, u_act,'--','LineWidth',lw)
xlabel('x direction')
ylabel('u')
subplot(1,3,3)
hold on
plot(xe, u_new(2:end-1).*Hu,'LineWidth',lw)
plot(xc, u_act_c.*H_act,'--','LineWidth',lw)
xlabel('x direction')
ylabel('u')

