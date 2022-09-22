%Based on Savant et al (2018)
%Sloshing example with no friction

%VARIABLES
%H:   total elevatation of water, H = eta + B
%eta: water elevation above the plane y = 0
%B:   depth of terrain below the plane y = 0
%u:   depth-averaged velocity in x-direction

addpath('~/Documents/MATLAB/SWE-FVM/swe-1d/include')

clear all
close all

%parameters
t_f = 24.0 * 60 * 60; %[s]
Lx  = 4.0 * 10^4;     %[m]
g   = 9.81;           %[m/s^2]
a   = 0.1;            %[m]
H0  = 10.0;           %[m]

%boundary and initial data
C = 0.0; %const of int for act soln
[~,g_u] = act_solution(0,0,a,C,g,H0,Lx);    %u(0,0), not BC for all times
g_eta = act_solution(Lx,0,a,C,g,H0,Lx); %eta(Lx,0), not BC for all times

%Tolerances
n_tol  = 1.0e-10;
wd_tol = 1.0e-10;
max_iter = 10;

%timestepping parameters
u_max = g_u + sqrt(g*H0);
CFL = 2.0;
theta = 0.5;

%number of cells in each direction
N = 2^6;

%step sizes
dx = Lx/N;
dt = CFL * dx / u_max;
t_incr = 500;
stride = ceil(t_incr/dt);

nu = 0.00*dx^2; %numerical viscosity

%cell centers
xc = ((1:N)' - 0.5) * dx;

%edges 
xe = (0:N)' * dx;

%initialize total water column, bathymetry, and surface elevation
eta = act_solution(xc,0,a,C,g,H0,Lx);
B   = H0*ones(N,1);
H   = eta + B;
      
%axis
eta_axis = [-2*dx Lx+2*dx -1.5*a 1.5*a];
u_axis   = [-2*dx Lx+2*dx -0.5*u_max*a 1.5*u_max*a];
q_axis   = [-2*dx Lx+2*dx -0.5*u_max*H0*a 1.5*u_max*H0*a];

%initialize bathymetry
B = [B(1); B; B(end)];

%initialize variables
[~,u] = act_solution(xe,0,a,C,g,H0,Lx);

%create ghost cell arrays
eta_new = [eta(1); eta; eta(end)];
H_new = eta_new + B;
u_new = [g_u; u; u(end)];

%intiialize variables for previous time steps
eta_old = eta_new;
H_old = H_new;
u_old = u_new;

u_old_old = u_old;

gamma = zeros(N+1,1);
Hu = calculate_edge_elevations(u_old,H_old);

%to track oscillations at specified x location
x0 = Lx/2;
x_idx = ceil(x0/dx);
num_ts = ceil(t_f/dt);
pos    = zeros(num_ts,1);

t_idx = 0;
t = 0;
while t <= t_f
              
    %tracking oscillations
    pos(t_idx+1) = eta_old(x_idx+1);
    
    %set ghost cells
    [g_eta,~] = act_solution(Lx,t,a,C,g,H0,Lx);
    eta_old = set_eta_gc(eta_old,g_eta);
    u_old = set_u_gc(u_old);
        
    %plot 
    plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct u at non-grid points for semi-Lagrangian approach
    u_star = construct_u_star(u_old,u_old_old,dx,dt);
    
    %gradient of eta at t_n
    d_eta_n_dx = (eta_old(2:end) - eta_old(1:end-1))/dx;
    
    %Laplacian of u,v at t_n
    Lap_u = Lap_1d(u_star,dx);

    %Casulli's "G" operator
    Gu = construct_Gu(u_star,Lap_u,d_eta_n_dx,g,nu,dt,theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_old = eta_old + B;
    Hu = calculate_edge_elevations(u_star,H_old);
      
    %wet vector
    w = find(Hu > wd_tol);  
    
    %friction
    gamma = zeros(N+1,1);    
    
    %prefactor used for T and velocity update
    PHI = zeros(N+1,1);
    PHI(w) = Hu(w) ./ (Hu(w) + dt*gamma(w));
    
    %construct V old (linear)
    V_old = H_old * dx;
    V_old(2:end-1) = total_vol(xc,eta_old(2:end-1),dx,wd_tol);
    
    %construct fluxes
    [flux_x_ex,flux_x_im] = construct_fluxes(Hu,u_old,Gu,PHI);
    [flux_x_ex,flux_x_im] = set_flux_bc(flux_x_ex,flux_x_im,Hu,u_old,Gu,PHI,g_u);
        
    %update volume
    V_new = V_old(2:end-1) - dt*flux_x_ex;
    
    %reconstruct total elevation H
    H_new(2:end-1) = V_new / dx;
    
    %calculate surface elevation eta
    eta_new(2:end-1) = H_new(2:end-1) - B(2:end-1);
    
    %nonlinear
    T = construct_T(Hu,eta_old,theta,dt,dx,g,PHI);
    b = construct_b(V_old,flux_x_ex,flux_x_im,dt,theta);
    [eta_update, r,iter] = nonlinear_solve(xc,dx,eta_old,T,b,n_tol,wd_tol,max_iter);
    eta_new(2:end-1) = eta_update;
    
    eta_new = set_eta_gc(eta_new,g_eta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_eta_n1_dx = (eta_new(2:end) - eta_new(1:end-1))/dx;
    
    %update u
    RHS = Gu - dt*g*theta*d_eta_n1_dx;
    u_new(2:end-1) = PHI.*RHS;
     
    %set u bc's
    [~,g_u] = act_solution(0,t,a,C,g,H0,Lx);
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
[eta_act,u_act_c] = act_solution(xc,t,a,C,g,H0,Lx);
H_act = eta_act + H0;
q_act = H_act .* u_act_c;
[~,u_act] = act_solution(xe,t,a,C,g,H0,Lx);

%error in SS soln
L_2_e  = [norm(eta_act - eta_new(2:end-1))     norm(u_act - u_new(2:end-1))] * sqrt(dx)
L_oo_e = [max(abs(eta_act - eta_new(2:end-1))) max(abs(u_act - u_new(2:end-1)))]


plot_state(xc,xe,t,0,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
subplot(1,3,1)
hold on
plot(xc, eta_act)
plot(xc,g_eta*ones(N,1),'--')
plot(xc,-g_eta*ones(N,1),'--')
plot(xc, -B(2:end-1))
plot(x0*ones(N,1),H0*(xc/Lx-0.5))
subplot(1,3,2)
hold on
plot(xe, u_act)
subplot(1,3,3)
hold on
plot(xc, q_act)

%plot oscillations at x0
figure(2)
hold on
ts = 0:dt:t_f;
pos_act = act_solution(xc(x_idx),ts,a,C,g,H0,Lx);
plot(ts,pos)
plot(ts,pos_act)
plot(ts,max(pos)*ones(num_ts,1))
ylabel('\eta(x_0,t)')
xlabel('t')
axis([0 t_f -a/2 a/2])
