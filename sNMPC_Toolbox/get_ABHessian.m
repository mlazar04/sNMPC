function [sys_new] = get_ABHessian(sys)
% Computes Jacobian and Hessian matrices for Taylor approxiamtions
sys_new=sys;
% Linearization around (0,0)
x0=zeros(sys.n,1);
u0=zeros(sys.m,1);

sys_new.A=double(subs(jacobian(sys.xdot,sys.x),[sys.x;sys.u],[x0;u0]));
sys_new.B=double(subs(jacobian(sys.xdot,sys.u),[sys.x;sys.u],[x0;u0]));

% Hessian to be used with quasi-second order approximation
Hessian=cell([1 sys.n]);
for i=1:sys.n
    sys_new.Hessian{1,i}=double(subs(hessian(sys.xdot(i),[sys.x;sys.u]),[sys.x;sys.u],[x0;u0]));
end
end