clear 
close all
clc
format compact 
format long

%% System equations in input affine form xdot=fx+fu*u

% cart spring-damper system
sys.Ts=0.4;
syms x1 x2 u
sys.x=[x1;x2];
sys.u=u;

rho0=0.33;
Mc=1;
hd=1.1;

sys.fx=[x1+sys.Ts*x2;-sys.Ts*rho0/Mc*exp(-x1)*x1+(1-sys.Ts*hd/Mc)*x2];
sys.fu=[0;sys.Ts/Mc];
sys.xdot=sys.fx+sys.fu*sys.u;

%% Constraints
sys.n=size(sys.x,1);
sys.m=size(sys.u,1); 

% Definition of state and input constraints
sys.x_low=[];
sys.x_high=[];
for i=1:2:sys.n-1
    sys.x_low=[sys.x_low;[-2.65;-5]];
    sys.x_high=[sys.x_high;[2.65;5]];
end
sys.u_low=-4.5*ones(sys.m,1);
sys.u_high=4.5*ones(sys.m,1);

%% Tuning parameters, optimizer settings 

% Number of terminal sets
p.M=1;
% Horizon
p.N=5;
% Cost matrices
p.Q=blkdiag(1,1)*0.01;
p.R=blkdiag(10);

% Display M N Q R
M = p.M
N = p.N
Q = p.Q
R = p.R

% Relaxation parameter in the nonlinear inequality
p.kappaj=0.5*ones(1,p.M);

% Scaling of the state-space to define LDI region
p.vscale=0.9;

% Parameter to enforce variation between terminal sets
p.rho=1.005;

% Number of random restarts for a possibly increased nonlinear optimization accuracy
p.gridpoints=30;

% Maximum scaling during bisection
p.a_bar_max=1e3;
% Tolerance during bisection
p.tol=1e-3;

% Nonlinear optimizer settings
opt_NL=sdpsettings('solver','ipopt','verbose',0);
opt_NL.usex0=1;

% Linear optimizer settings (change solver to 'sedumi' if MOSEK is not installed)
opt_L=sdpsettings('solver','mosek','verbose',0);

% MPT YSet settings (used only for plotting) (change solver to 'sedumi' if MOSEK is not installed)
ysetoptions=sdpsettings('solver','mosek','verbose',0);

% Type of the approxiamtion and contorl law
fprintf('\n1: Linear 2: Quasi-second order 3: Linear+NLcontrol 4: Quasi-second order+NLcontrol\n');
choice='Choose approximation order and terminal control law:';
Mode=input(choice);

%% Linearization and computation of the Hessian for LDI approximation
[sys]=get_ABHessian(sys);

%% Solve the LMIs 
[P, K, alpha, E1, VOL1, XUset, Xset_scaled]=solve_LMIs(sys,p,Mode,opt_L);

%% Plot resulting terminal sets before bisection
figure(); hold on;
xlabel('x1') 
ylabel('x2') 
plot_ellipsoidal_sets(sys, p, E1, VOL1, XUset, Xset_scaled);
hold off

%% Nonlinear optimization and bisection
[alphascale, E2, VOL2]=solve_nlp_bisection(sys,p,P,K,alpha,Mode,opt_NL);
max(cell2mat(VOL2)) % Print volume of the largest set

%% Plot resulting terminal sets after bisection
figure(); hold on
xlabel('x1') 
ylabel('x2') 
plot_ellipsoidal_sets(sys, p, E2, VOL2, XUset, Xset_scaled);
hold off

%% Simulation in CasADi
fprintf('\nWould you like to simulate the system or obtain DOA plot?\n' );
fprintf('CasADi is required (must be added to path)!\n')
choice = 'y/n:';

if input(choice,"s") == 'y'
    % import CasADi
    import casadi.*

    % Define states,inputs and system equations in CasADi format
    s.Ts=sys.Ts;
    s.x1=SX.sym('x1');
    s.x2=SX.sym('x2');
    s.u=SX.sym('u');

    s.x = [s.x1; s.x2];

    s.fx=[s.x1+s.Ts*s.x2;-s.Ts*rho0/Mc*exp(-s.x1)*s.x1+(1-s.Ts*hd/Mc)*s.x2];
    s.fu=[0;s.Ts/Mc];
    s.xdot=s.fx+s.fu*s.u;

    s.x_low=sys.x_low;
    s.x_high=sys.x_high;
    s.u_low=sys.u_low;

    s.u_high=sys.u_high;

    %% Plot domian of attraction
    x_axis = -2.65:0.15:2.65;
    y_axis = -5:0.25:5;

    fprintf('\nWould you like to obtain a domain of attraction representation?\n' );
    fprintf('This takes a few minutes! Only available in 2D!\n')
    choice = 'y/n:';

    if input(choice,"s") == 'y'
        [feasible_points, infeasible_points]=plot_DOA(s,sys,p,P,alpha,alphascale,x_axis,y_axis,E2, VOL2, XUset, Xset_scaled);
    end
    %% Find initial terminal set for an initial condition
    x0=[-2;4];
    [feasible,init_index]=find_init_set(s,p,P,alpha,alphascale,x0);        

    %% Simulate the system
    sim_tim = 20; % Maximum simulation time

    [traj,t,ss_error,cost]=casadi_simulation(s,p,P,K,alpha,alphascale,x0,init_index,sim_tim);
    length(traj)

    %% Plot state trajectories and NMPC cost evolution

    % State trajectory
    figure()
    plot_ellipsoidal_sets(sys, p, E2, VOL2, XUset, Xset_scaled);
    plot(x0(1),x0(2),'o','MarkerFaceColor','r')
    plot(traj(1,:),traj(2,:),'r','LineWidth',1);
    title('State trajectory and terminal sets')

    % Evolution of the NMPC cost function
    figure()
    hold on
    plot(0:length(cost)-1,cost,'r','LineWidth',1);
    title('Evolution of the NMPC cost function')
    xlabel('iteration') 
    ylabel('value of the cost function') 

    % Evolution of system states 
    figure()
    hold on; 
    plot(0:length(traj)-1,traj(1,:),'-','Color','b','LineWidth',1); 
    plot(0:length(traj)-1,traj(2,:),'--','Color','b','LineWidth',1);
    title('Evolution of system states')
    legend('x1','x2');
end