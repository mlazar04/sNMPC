%% Randomly select points inside the hyperellipse with a unfiorm distribution
function [z_fa] = unif_sample(P1,gammascale,gridpoints)

E=ellipsoid(gammascale*inv(P1));
S=inv(P1);
Gamma_Threshold=gammascale*10;
z_hat=0*ones(length(P1),1);
m_FA=gridpoints;

nz=length(S);
X_Cnz=randn(nz,m_FA);
X_Cnz=X_Cnz./kron(ones(nz,1),sqrt(sum(X_Cnz.^2))); % Points uniformly distributed on hypersphere
R=ones(nz,1)*(rand(1,m_FA).^(1/nz )); % Points with pdf nzrˆ(nz1); 0<r<1
unif_sph=R.*X_Cnz; % m FA points in the hypersphere
T=chol(S); % Cholesky factorization of S => S=T’T
unif_ell =T'*unif_sph ; % Hypersphere to hyperellipsoid mapping
z_fa=(unif_ell*sqrt(Gamma_Threshold)+(z_hat*ones(1,m_FA))); % Translation around gate cente
%plot(z_fa(1,:),z_fa(2,:),'*');hold off;

end