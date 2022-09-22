clear 
close all
clc
format compact 
format long

%% System equations in input affine form xdot=fx+fu*u

% oscillator
sys.Ts=0.02;
syms x1 x2 x3 x4 u
sys.x=[x1;x2;x3;x4];
sys.u=u;

Jl=1; Jr=1; M1=0.1; l1=1; g=9.8; k=100;  

sys.fx(1:2,1)=[x1+sys.Ts*x2;x2+sys.Ts*(-M1*g*l1/Jl*sin(x1)-k/Jl*(x1-x3))];
sys.fx(3:4,1)=[x3+sys.Ts*x4;x4+sys.Ts*(k/Jr*(x1-x3))];
sys.fu=sys.Ts*[0; 0; 0; 1/Jr];
sys.xdot=sys.fx+sys.fu*sys.u;

%% Constraints
sys.n=size(sys.x,1);
sys.m=size(sys.u,1); 

% Definition of state and input constraints
sys.x_low=[-pi; -pi; -10; -10];
sys.x_high=[pi; pi; 10; 10];
sys.u_low=-2*ones(sys.m,1);
sys.u_high=2*ones(sys.m,1);

%% Tuning parameters, optimizer settings 

% Number of terminal sets
p.M=1;
% Horizon
p.N=5;
% Cost matrices
p.Q=blkdiag(10,1,10,1)*0.1;
p.R=0.001;

% Display M N Q R
M = p.M
N = p.N
Q = p.Q
R = p.R

% Relaxation parameter in the nonlinear inequality
p.kappaj=0.0075*ones(p.M);

% Scaling of the state-space to define LDI region
p.vscale=1;

% Parameter to enforce variation between terminal sets
p.rho=1.001;

% Number of random restarts for a possibly increased nonlinear optimization accuracy
p.gridpoints=20;

% Maximum scaling during bisection
p.a_bar_max=1e3;
% Tolerance during bisection
p.tol=1e-3;

% Nonlinear optimizer settings (change solver to 'fmincon' if IPOPT is not installed)
opt_NL=sdpsettings('solver','ipopt','verbose',0);
opt_NL.usex0=1;

% Linear optimizer settings (change solver to 'sedumi' if MOSEK is not installed)
opt_L=sdpsettings('solver','mosek','verbose',0);

% MPT YSet settings (used only for plotting) (change solver to 'sedumi' if MOSEK is not installed)
ysetoptions=sdpsettings('solver','mosek','verbose',0);

% Type of the approxiamtion and contorl law
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
fprintf('\nWould you like to simulate the system?\n' );
fprintf('CasADi is required (must be added to path)!\n')
choice = 'y/n:';

if input(choice,"s") == 'y'

    % import CasADi
    import casadi.*

    % Define states,inputs and system equations in CasADi format
    s.Ts=sys.Ts;
    s.x1=SX.sym('x1');
    s.x2=SX.sym('x2');
    s.x3=SX.sym('x3');
    s.x4=SX.sym('x4');
    s.u=SX.sym('u');

    s.x = [s.x1; s.x2; s.x3; s.x4];

    s.fx=[s.x1+s.Ts*s.x2;s.x2+s.Ts*(-M1*g*l1/Jl*sin(s.x1)-k/Jl*(s.x1-s.x3));
          s.x3+s.Ts*s.x4;s.x4+s.Ts*(k/Jr*(s.x1-s.x3))];
    s.fu=s.Ts*[0; 0; 0; 1/Jr];
    s.xdot=s.fx+s.fu*s.u;

    s.x_low=sys.x_low;
    s.x_high=sys.x_high;
    s.u_low=sys.u_low;
    s.u_high=sys.u_high;


    %% Find initial terminal set for an initial condition
    x0=[-1;1;-1;1];
    [feasible,init_index]=find_init_set(s,p,P,alpha,alphascale,x0)

    %% Simulate the system
    sim_tim = 20; % Maximum simulation time

    [traj,t,ss_error,cost]=casadi_simulation(s,p,P,K,alpha,alphascale,x0,init_index,sim_tim);
    length(traj)

    %% Plot state trajectories and NMPC cost evolution

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
    plot(0:length(traj)-1,traj(3,:),'-','Color','r','LineWidth',1);
    plot(0:length(traj)-1,traj(4,:),'--','Color','r','LineWidth',1);
    title('Evolution of system states')
    legend('x1','x2','x3','x4');
end