%% Construct and solve the LMIs
function [P, K, alpha, E1, VOL1, XUset, Xset] = solve_LMIs(sys, p, Mode, opt_L)

% Check if the Ellipsoidal Toolbox can be used
try
    testET=ellipsoid(eye(2));
    ET=1;
catch
    ET=0;
end

% Set mode dependent parameters
if Mode==2||Mode==4
    ldilmi=1;
else
    ldilmi=0;
    vscale=1;
end

% Assign system parametrs to local variables
n=sys.n; m=sys.m; A=sys.A; B=sys.B; Hessian=sys.Hessian;
u_high=sys.u_high; u_low=sys.u_low; x_high=sys.x_high; x_low=sys.x_low;
M=p.M; Q=p.Q; R=p.R;
vscale=p.vscale;
rho=p.rho;
a_bar_max=p.a_bar_max;

% Use a constant kappa 
kappa=cell([1 M]);
for i=1:M
    kappa{1,i}=p.kappaj(M);
end

% Set input and state admissable sets
Uset=Polyhedron([diag(1./u_high);diag(1./u_low)],[ones(size(u_high,1),1);ones(size(u_high,1),1)]);
Uset=Polyhedron((Uset.A./repmat(Uset.b, 1, Uset.Dim)),ones(size(Uset.b)));
Xset=Polyhedron([diag(1./x_high);diag(1./x_low)],[ones(size(x_high,1),1);ones(size(x_high,1),1)]);
Xset=Polyhedron( (Xset.A./repmat(Xset.b, 1, Xset.Dim)),ones(size(Xset.b)));

% Initiate empty cell arrays
O=cell([1 M]);
Y=cell([1 M]);
P=cell([1 M]);
K=cell([1 M]);
XUset=cell([1 M]);
E1=cell([1 M]);
VOL1=cell([1 M]);

% Try different scalings of the Xset for the quasi-second order approximation case
for scalingv=(vscale:-vscale/30:vscale/10)*ldilmi+(1-ldilmi)
    xldi=(scalingv*Xset.V)';
    Xset_scaled=Polyhedron('V',xldi');
    Xset_scaled=Polyhedron((Xset_scaled.A./repmat(Xset_scaled.b, 1, Xset_scaled.Dim)),ones(size(Xset_scaled.b)));
    
    L=[];
    yalmip('clear')
    sdpvar a_bar

%   Define variables
    for i=1:M
        O{1,i}=sdpvar(n,n);
        Y{1,i}=sdpvar(m,n);
    end
    
%   Add constraints to the linear optimization
%   Lyapunov function 1 to M-1
    for i=1:M-1
%   Vertices of the scaled polytope for LDI case
        for v=1:size(xldi,2)*ldilmi+(1-ldilmi)
            Gxx=[];
            Gxu=[];
%           Build the quasi-second order approximation
            for it=1:n
                Gxx=[Gxx;0.5*double(xldi(:,v)'*Hessian{1,it}(1:n,1:n))];
                Gxu=[Gxu;0.5*double(xldi(:,v)'*Hessian{1,it}(1:n,n+1:n+m))];
            end
            
            Abar=A+Gxx*ldilmi;
            Bbar=B+Gxu*ldilmi;
            
            L1=[[(1-kappa{1,i})*O{1,i} (Abar*O{1,i}+Bbar*Y{1,i})' O{1,i} Y{1,i}';
                (Abar*O{1,i}+Bbar*Y{1,i}) O{1,i+1} zeros(n,n) zeros(n,m);
                O{1,i} zeros(n,n) a_bar*inv(Q) zeros(n,m);
                Y{1,i} zeros(m,n) zeros(m,n) a_bar*inv(R)]>=1e-6*eye(3*n+m)]:'Lyapunov function';
            L=[L L1];
        end
        L2=[O{1,i}>=1e-6*eye(n)]:'Positive definite';
        L=[L L2];    
    end
    
%   Lyapunov function M comes back to 1
    for v=1:size(xldi,2)*ldilmi+(1-ldilmi)
        Gxx=[];
        Gxu=[];
        for it=1:n
            Gxx=[Gxx;0.5*double(xldi(:,v)'*Hessian{1,it}(1:n,1:n))];
            Gxu=[Gxu;0.5*double(xldi(:,v)'*Hessian{1,it}(1:n,n+1:n+m))];
        end
            
        Abar=A+Gxx*ldilmi;
        Bbar=B+Gxu*ldilmi;
        L1=[[(1-kappa{1,M})*O{1,M} (Abar*O{1,M}+Bbar*Y{1,M})' O{1,M} Y{1,M}';
            (Abar*O{1,M}+Bbar*Y{1,M}) O{1,1} zeros(n,n) zeros(n,m);
            O{1,M} zeros(n,n) a_bar*inv(Q) zeros(n,m);
            Y{1,M} zeros(m,n) zeros(m,n) a_bar*inv(R)]>=1e-6]:'Lyapunov function M';
        L=[L L1];
    end
    L2=[O{1,M}>=1e-6]:'Positive definite M';
    L=[L L2];
    
%   Force variation within the set sequence
    for i=1:M-1
        L3=[Y{1,i+1} >= rho*Y{1,i}]:'scaling1';
        L=[L L3];
    end
    
 %  State and input constraints
    for i=1:M
        for j=1:size(Xset_scaled.A,1)
            L4=[Xset_scaled.A(j,:)*O{1,i}*Xset_scaled.A(j,:)'<=1]:'LDI constraint';
            L=[L L4];
        end
        for j=1:size(Uset.A,1)
            L5=[[1 (Uset.A(j,:))*Y{1,i};Y{1,i}'*(Uset.A(j,:))' O{1,i}]>=0]:'input constraint' ;
            L=[L L5];
        end
    end
    
%   Bound a_bar
    Le=[1e-6<=a_bar a_bar<=a_bar_max];
    L=[L Le];

%   Solve the LMIs while minimizing over -log(det(O)
    o=optimize(L,-log(det(O{1,max(floor(M/2),1)})),opt_L)
    %o=optimize(L,-log(det(O{1,1})),opt_L)

    if o.problem==0
        for i=1:M
%           Compute P and K matrices
            alpha=1/(value(a_bar));
            P{1,i}=inv(value(O{1,i}));
            K{1,i}=value(Y{1,i})*inv(value(O{1,i}));
            
%           Construct the combined state and input constrained set
            Kconstruct=Uset.A*K{1,i};
            XUset{1,i}=Polyhedron([Xset.A;Kconstruct],[Xset.b;Uset.b]);
            XUset{1,i}=Polyhedron((XUset{1,i}.A./repmat(XUset{1,i}.b, 1, XUset{1,i}.Dim)),ones(size(XUset{1,i}.b)));
            
%           Construct ellipsoids and compute their volumes
            if ET==1
                E1{1,i}=ellipsoid(inv((P{1,i})));
                VOL1{1,i}=volume(E1{1,i});
            else
                x_temp=sdpvar(n,1);
                E1{1,i}=YSet(x_temp,x_temp'*P{1,i}*x_temp<=1,opt_L);
                VOL1{1,i}=sqrt(det(inv(P{1,i})));
            end
        end
        
%       If problem is feasible no need to further sclae down the polytope
        scalingv
        break
        
    else
%       Continue with further scaling down the polytope
        if Mode==2||Mode==4
            check(L)
            disp(['Infeasible for scalingv=', num2str(scalingv),'.Further scaling down continues...'])
        else
            check(L)
            error('Infeasible')
            break;
        end
    end   
end

P{1,M+1}=P{1,1};
K{1,M+1}=K{1,1};

end
