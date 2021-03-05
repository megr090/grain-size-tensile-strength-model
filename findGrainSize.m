function [grainsize,tau,A] = findGrainSize(T,theta,strainrate,R,mu,D,M0,k0,gamma,p,c,n)
% This function computes grain size from temperature

[Qg,Qc,Qm] = defineActivationEnergies(T);

A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
A = zeros(size(T));
for i=1:length(T)
    A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
end
tau = A.^(-1/n).*strainrate.^(1/n); % Pa

% compute grain size
grainsize = zeros(1,length(T));
grainsize = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T)))./(4.*mu.^2.*tau.*(1-theta).*strainrate)).^(1/(1+p));
grainsize =  grainsize.*1e3;

end
