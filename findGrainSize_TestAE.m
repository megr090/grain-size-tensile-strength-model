function [grainsize,viscosity,time_to_eq,tau,A] = findGrainSize_TestAE(T,f,strainrate,R,mu,D,M0,k0,gamma,p,c,n,Qgm,Qgp,tc)

[Qg,Qc,Qm] = defineActivationEnergies_Test(T,Qgm,Qgp,tc);

A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
A = zeros(size(T));
for i=1:length(T)
    A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
end
tau = A.^(-1/n).*strainrate.^(1/n); % Pa

% compute grain size
grainsize = zeros(1,length(T));
%grainsize_temp(i) = ((k0.*exp(-Qg(i)./(R.*T(i))).*p.^(-1).*cons1.*gamma)./(tau(i).*(1-f).*epsilon_disl)).^(1/(1+p));
grainsize = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T)))./(4.*mu.^2.*tau.*(1-f).*strainrate)).^(1/(1+p));
grainsize =  grainsize.*1e3;

viscosity = zeros(1,length(T))+0.5.*(A.*strainrate.^(n-1)).^(-1/n); % Pa s

% compute equilibrium time
time_to_eq = (p.*(1/0.05).*strainrate).^(-1); % in s
time_to_eq = time_to_eq./(365*24*60*60);

end