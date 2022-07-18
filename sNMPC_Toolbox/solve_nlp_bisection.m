%% Checks the nonlinear inequality and scalaes down the sets by bisecito
function [gammascale, E, VOL2] = solve_nlp_bisection(sys, p, P, K, alpha, Mode, opt_NL)

% Assign parameters to local variables
A=sys.A; B=sys.B; Hessian=sys.Hessian; n=sys.n; m=sys.m;
x=sys.x; u=sys.u; fx=sys.fx; fu=sys.fu; xdot=sys.xdot;
x_high=sys.x_high; x_low=sys.x_low;
M=p.M; Q=p.Q; R=p.R; tol=p.tol; gridpoints=p.gridpoints; kappaj=p.kappaj;
ALPHA=alpha;

% Initialize cell arrays
uT=cell([1 M]);
r1=cell([1 M]);
r2=cell([1 M]);
r2_fun=cell([1 M]);
r_fun=cell([1 M]);
Abar=cell([1 M]);
Bbar=cell([1 M]);
Abar_fun=cell([1 M]);
Bbar_fun=cell([1 M]);
VOL2=cell([1 M]);

% Set mode dependent parameters
if Mode==2||Mode==4
    ldilmi=1;
else
    ldilmi=0;
    vscale=1;
end

for k=1:M
%   Define terminal controll law
%   1: Linear 2: Quasi-second 3: Linear+NLcontrol 4: Quasi-second+NLcontrol
    if Mode==1||Mode==2
%       Linear state feedback
        uT{1,k}=K{1,k}*x;
    else
%       Nonlinear feedback
        uT{1,k}=-(fu'*P{1,k}*fu+R*ALPHA)\(fu'*P{1,k}*fx);
        %uT{1,k}=simplifyFraction(uT{1,k});
    end
    
%   Define r2(x) as a function of x
    r2{1,k}=uT{1,k}-K{1,k}*x;
    r2_fun{1,k}=matlabFunction(r2{1,k},'Vars',{x});

    
%   Quasi-second order approximation
    Gxx=[];
    Gxu=[];
    for it=1:n
        Gxx=[Gxx;0.5*(x'*Hessian{1,it}(1:n,1:n))];
        Gxu=[Gxu;0.5*(x'*Hessian{1,it}(1:n,n+1:n+m))];
    end
    
    % Construct approxiamtion matrices
    Abar{1,k}=A+Gxx*ldilmi;
    Abar_fun{1,k}=matlabFunction(Abar{1,k},'Vars',{x});
    
    Bbar{1,k}=B+Gxu*ldilmi;
    Bbar_fun{1,k}=matlabFunction(Bbar{1,k},'Vars',{x});
    
    % Define r1(x) 
    r1{1,k}=subs(xdot,u,K{1,k}*x)-(Abar{1,k}*x+Bbar{1,k}*K{1,k}*x);
    
    % Define the residual r(x) as a function of x
    r_fun{1,k}=matlabFunction(r1{1,k}+fu*r2{1,k},'Vars',{x});
end

%% Nonlinear optimization and bisection
gammascale=1;
iter=0;
x_temp=sdpvar(n,1);

fprintf('\nstart nonlinear inequality check\n');
while 1
    iter=iter+1
    
    yalmip('clear')
    clear Objective Constraints
    
    x_temp=sdpvar(n,1);
    
    P1=P{1,iter};
    P2=P{1,iter+1};
    kappai=kappaj(iter);
    r=r_fun{1,iter}(x_temp);
    r22=r2_fun{1,iter}(x_temp);
    Ki=K{1,iter};
    ABAR=Abar_fun{1,iter}(x_temp);
    BBAR=Bbar_fun{1,iter}(x_temp);
    
%   Create random points within the hyperellipsoid
    xrand = unif_sample(P1,gammascale,gridpoints);
    xrand(:,1) = zeros(n,1);
    
%   fopt is assigning random initial conditions for the nonlinear optimization
    for fopt=1:gridpoints
        % Define NLP
        assign(x_temp,xrand(:,fopt));
        Objective=-(2*x_temp'*(ABAR+BBAR*Ki)'*P2*r+r'*P2*r+2*x_temp'*Ki'*R*ALPHA*r22+r22'*R*ALPHA*r22-kappai*x_temp'*P1*x_temp);
        Constraints=[x_temp'*P1*x_temp<=gammascale x_low<=x_temp x_temp<=x_high];
        
        o=optimize(Constraints,Objective,opt_NL);
        if o.problem==1
            error('');
        end
        Rmax=-1*value(Objective);
        
%       If any point found breaking the upperbound start bisection
        if (Rmax>0)
            break
        end
    end
    
    if ~(Rmax>0)
%       Inequality satisifed no need for further bisection
        if iter==max(M,1)
            break
        end
    else
        fprintf('nonlinear inequality not satisfied: bisection\n');
        gammascale_upper=gammascale;
        gammascale_lower=eps;
        
%       Bisection        
        while (gammascale_upper-gammascale_lower)>(tol)
            gammascale_test = (gammascale_upper+gammascale_lower)/2;
            
            yalmip('clear')
            clear Objective Constraints
            x_temp=sdpvar(n,1);
            
            r=r_fun{1,iter}(x_temp);
            r22=r2_fun{1,iter}(x_temp);
            ABAR=Abar_fun{1,iter}(x_temp);
            BBAR=Bbar_fun{1,iter}(x_temp);
            
            %[xrand] = unif_sample(P1,gammascale_test,gridpoints);
            for fopt=1:gridpoints
                % Define NLP
                assign(x_temp,xrand(:,fopt));
                Objective=-(2*x_temp'*(ABAR+BBAR*Ki)'*P2*r+r'*P2*r+2*x_temp'*Ki'*R*ALPHA*r22+r22'*R*ALPHA*r22-kappai*x_temp'*P1*x_temp);
                Constraints=[x_temp'*P1*x_temp<=gammascale_test x_low<=x_temp x_temp<=x_high];
                
                o=optimize(Constraints,Objective,opt_NL);
                if o.problem==1
                    error('');
                end
                
                Rmax=-1*value(Objective);
%           If any point found breaking the upperbound back to bisection
                if (Rmax>0)
                    break
                end
            end
            
%           Update bisection bounds
            if Rmax>0
                gammascale_upper=gammascale_test;
            else
                gammascale_lower=gammascale_test;
            end
        end
        if gammascale_test==eps
            fprintf('\nproblem?\n');
            error('scaling error');
        end
        
%       Update gammascale
        gammascale=gammascale_test
        if iter==max(M,1)
            break
        end
        iter
    end
end



%% Update the ellipsoids and volumes
try
    testET=ellipsoid(eye(2));
    ET=1;
catch
    ET=0;
end

for i=1:M
    if ET==1
        E{1,i}=ellipsoid(inv((P{1,i}/gammascale)));
        VOL2{1,i}=volume(E{1,i});
    else
        x_temp=sdpvar(n,1);
        E{1,i}=YSet(x_temp,x_temp'*P{1,i}*x_temp<=gammascale,ysetoptions);
        VOL2{1,i}=sqrt(det(inv(P{1,i}/gammascale)));
    end
end
end