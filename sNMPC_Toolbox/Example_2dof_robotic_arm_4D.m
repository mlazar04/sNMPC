clear 
close all
clc
format long

%% System equations in input affine form xdot=fx+fu*u

% 2DOF robotic arm
sys.Ts=0.01;
syms x1 x2 x3 x4 u1 u2
sys.x=[x1;x2;x3;x4];
sys.u=[u1;u2];

b1 = 0.0715; b2 = 0.0058; b3 = 0.0114; b4 = 0.3264; b5 = 0.3957;
b6 = 0.6254; b7 = 0.0749; b8 = 0.0705; b9 = 1.1261; 

M1 = [b1+b2+2*b3*cos(x2) b2+b3*cos(x2); b7-b8*cos(x2) b7];
c = [-b3*x4^2*sin(x2)-2*b3*x3*x4*sin(x2)+b6*x3; b3*x3^2*sin(x2)+b9*x4];
g = [-b4*sin(x1+x2)-b5*sin(x1); -b4*sin(x1+x2)];


fx(1:2,1)=[x1+sys.Ts*x3;x2+sys.Ts*x4];
fx(3:4,1)=[x3;x4]-sys.Ts*(inv(M1)*(c+g));
sys.fu=[0 0; 0 0; sys.Ts*inv(M1)];
sys.fx=simplify(fx);
sys.xdot=sys.fx+sys.fu*sys.u;

%% Constraints
sys.n=size(sys.x,1);
sys.m=size(sys.u,1); 

% Definition of state and input constraints
sys.x_low=[-pi/2; -pi/2; -10; -10];
sys.x_high=[pi/2; pi/2; 10; 10];
sys.u_low=-1*ones(sys.m,1);
sys.u_high=1*ones(sys.m,1);

%% Tuning parameters, optimizer settings 

% Number of terminal sets
p.M=1;
% Horizon
p.N=5;
% Cost matrices
p.Q=blkdiag(10,10,1,1);
p.R=blkdiag(1,1);

% Display M N Q R
M = p.M
N = p.N
Q = p.Q
R = p.R

% Relaxation parameter in the nonlinear inequality
p.kappaj=0.005*ones(1,p.M);

% Scaling of the state-space to define LDI region
p.vscale=1;

% Parameter to enforce variation between terminal sets
p.rho=1.05;

% Number of random restarts for a possibly increased nonlinear optimization accuracy
p.gridpoints=2;

% Maximum scaling during bisection
p.a_bar_max=1e6;
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
fprintf('\n1: Linear 2: LDI 3: Linear+NLcontrol 4: LDI+NLcontrol\n');
choice='Choose scaling via nonlinear optimisation method:';
Mode=input(choice);

%% Linearization and computation of the Hessian for LDI approximation
[sys]=get_ABHessian(sys);

%% Solve the LMIs 
[P, K, alpha, E1, VOL1, XUset, Xset_scaled]=solve_LMIs(sys,p,Mode,opt_L);
max(cell2mat(VOL1))

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

% Add CasADi to your MATLAB path
import casadi.*

% Define states,inputs and system equations in CasADi format
s.Ts=sys.Ts;
s.x1=SX.sym('x1');
s.x2=SX.sym('x2');
s.x3=SX.sym('x3');
s.x4=SX.sym('x4');
s.u1=SX.sym('u1');
s.u2=SX.sym('u2');

s.x = [s.x1; s.x2; s.x3; s.x4];
s.u = [s.u1; s.u2];

M1 = [b1+b2+2*b3*cos(s.x2) b2+b3*cos(s.x2); b7-b8*cos(s.x2) b7];
c = [-b3*s.x4^2*sin(s.x2)-2*b3*s.x3*s.x4*sin(s.x2)+b6*s.x3; b3*s.x3^2*sin(s.x2)+b9*s.x4];
g = [-b4*sin(s.x1+s.x2)-b5*sin(s.x1); -b4*sin(s.x1+s.x2)];

s.fx=[s.x1+s.Ts*s.x3;s.x2+s.Ts*s.x4; [s.x3;s.x4]-s.Ts*inv(M1)*(c+g)];
s.fu=[0 0; 0 0; M1/s.Ts]
s.xdot=s.fx+s.fu*s.u;

s.x_low=sys.x_low;
s.x_high=sys.x_high;
s.u_low=sys.u_low;
s.u_high=sys.u_high;

%% Find initial terminal set for an initial condition
x0 = [1.5; -1; 3; -7]; % initial condition
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
plot(traj(1,:),'-','Color','b','LineWidth',1); 
plot(traj(2,:),'--','Color','b','LineWidth',1);
plot(traj(3,:),'--','Color','r','LineWidth',1);
plot(traj(4,:),'--','Color','r','LineWidth',1);
title('Evolution of system states')
legend('x1','x2','x3','x4');