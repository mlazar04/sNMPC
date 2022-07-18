%% Simulate the NMPC for the obtained terminal ingredients
function [trajectory,time,ss_error,cost,u_cl]=casadi_simulation(s,p,P,K,alpha,gammascale,x0,init_index,sim_tim)
%%
import casadi.*

x=s.x;u=s.u;xdot=s.xdot;Ts=s.Ts;
x_low=s.x_low;x_high=s.x_high;u_low=s.u_low;u_high=s.u_high;
M=p.M;N=p.N;Q=p.Q;R=p.R;
n=length(s.x); m=length(s.u);
A_bar=1/alpha;

f = Function('f',{x,u},{xdot}); % nonlinear mapping function f(x,u)
U = SX.sym('U',m,N); % decision variables (controls)states over the optimization problem
X = SX.sym('X',n,(N+1)); % states over the optimization problem
Param = SX.sym('Param',n);

% Compute solution symbolically
X(:,1) = Param(1:n); % initial state
for i = 1:N
    st = X(:,i);  in = U(:,i);
    st_next  = f(st,in);
    X(:,i+1) = st_next;
end
% Function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,Param},{X});

% compute objective function for each terminal cost and constraint

for i=1:M
    obj = 0; % Objective function
    g = [];  % constraints vector
    for k=1:N
        st = X(:,k);  in = U(:,k);
        obj = obj + st'*Q*st + in'*R*in; % calculate obj
        g = [g ; X(:,k)];   %state xy
    end
    objective{1,i} = obj+A_bar*X(:,N+1)'*P{1,i}*X(:,N+1);
    constraints{1,i}=[g; X(:,N+1)'*P{1,i}*X(:,N+1)];
end

% make the decision variables one column vector
OPT_variables = reshape(U,m*N,1);
nlp_prob = struct('f', objective, 'x', OPT_variables, 'g', constraints, 'p', Param);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.error_on_fail = false;
%opts.ipopt.acceptable_tol =1e-8;
%opts.ipopt.acceptable_obj_change_tol = 1e-6;

args = struct;
% state constraints
args.lbg = [repmat(x_low,N,1); 0];  
args.ubg = [repmat(x_high,N,1); gammascale]; 

% input constraints
args.lbx = repmat(u_low,N,1);
args.ubx = repmat(u_high,N,1); 


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,m);  % two control inputs

% Start MPC
mpciter = 0;
u_cl=[];
cost=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm(x0,2) > 1e-2 && mpciter < sim_tim / Ts)
    args.p = x0; % set the values of the parameters vector
    args.x0 = reshape(u0',m*N,1); % initial value of the optimization variables
    %tic
    mode = mod(init_index + mpciter,M);
    if mode == 0
        mode = M;
    end
    solver = nlpsol('solver','ipopt',nlp_prob(mode), opts);
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
 
    %toc
    cost = [cost; full(sol.f)];
    u = reshape(full(sol.x)',m,N);
    %u = reshape(full(sol.x)',m,N)';
    x_eval = ff(u,args.p);
    x_term =full(x_eval(:,end));
    
    u = u';
    u_cl = [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = shift(Ts, t0, x0, u, f, K{1,mode},x_term); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0;
    mpciter = mpciter + 1;
end
trajectory=xx;
time = toc(main_loop);
ss_error = norm(x0,2);

end
