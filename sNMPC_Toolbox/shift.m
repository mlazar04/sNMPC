function [t0, x0, u0] = shift(T, t0, x0, u, f, K, x)
st = x0;
con = u(1,:)';
st = f(st,con);
x0 = full(st);

t0 = t0 + T;

u0 = [u(2:size(u,1),:); (K*x)'];
end