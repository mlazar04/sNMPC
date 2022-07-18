%% Determine initial pair of terminal ingrdients
function [feasible,ind] = find_init_set(s, p, P, alpha, gammascale, x0)
import casadi.*

x=s.x;u=s.u;xdot=s.xdot;
x_low=s.x_low;x_high=s.x_high;u_low=s.u_low;u_high=s.u_high;
M=p.M;N=p.N;Q=p.Q;R=p.R;
n=length(s.x); m=length(s.u);
A_bar=1/alpha;

f = Function('f',{x,u},{xdot}); % Nonlinear mapping function f(x,u)
U = SX.sym('U',m,N); % Decision variables (controls)states over the optimization problem
X = SX.sym('X',n,(N+1)); % States over the optimization problem
Param = SX.sym('Param',n);
gamma = SX.sym('gamma',M,1);

% Compute solution symbolically
X(:,1) = Param(1:n); % Initial state
for i = 1:N
    st = X(:,i);  in = U(:,i);
    st_next  = f(st,in);
    X(:,i+1) = st_next;
end
% Function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,Param},{X});

% Compute objective function for each terminal cost and constraint


obj = 0; % Objective function
g = [];  % Constraints vector
for k=1:N
    st = X(:,k);  in = U(:,k);
    obj = obj + st'*Q*st + in'*R*in; % Calculate obj
    g = [g ; X(:,k)];  % state x
end

% Define  convex combination of the P matrices
convcombP=0;
for i=1:M
    convcombP=convcombP+gamma(i)*P{1,i};
end
objective = obj+A_bar*X(:,N+1)'*convcombP*X(:,N+1); % Add terminal cost
constraints = [g; X(:,N+1)'*convcombP*X(:,N+1); sum(gamma)]; % Add terminal constraint

% Make the decision variables one column vector
OPT_variables = reshape(U,m*N,1);
nlp_prob = struct('f', objective, 'x', [OPT_variables; gamma], 'g', constraints, 'p', Param);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.error_on_fail = true;
%opts.ipopt.acceptable_tol =1e-8;
%opts.ipopt.acceptable_obj_change_tol = 1e-6;

args = struct;
% state constraints
args.lbg = [repmat(x_low,N,1); 0; 1];  
args.ubg = [repmat(x_high,N,1); gammascale; 1]; 

% input constraints
args.lbx = [repmat(u_low,N,1); zeros(M,1)];
args.ubx = [repmat(u_high,N,1); ones(M,1)];

%-------------------------------------------
u0 = zeros(N,m);  % m control inputs

args.p = x0; % set the values of the parameters vector
args.x0 = [reshape(u0',m*N,1); ones(M,1)]; % initial value of the optimization variables
solver = nlpsol('solver','ipopt',nlp_prob, opts);
feasible=false;
term_cost=[];
ind=0;

try
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(1:end-M))',m,N);
    xx = ff(u,x0);
    term_state=full(xx(:,end));

    for i=1:M
       term_cost=[term_cost term_state'*P{1,i}*term_state];
    end

    feasible=true;
    [min_cost,ind]=min(term_cost);

catch
    feasible=false;
    %fprintf('WARNING: Initial condition is not feasible!!!\n')
end