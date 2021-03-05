function [A,gsize,temp] = findPropertiesofShearMargin(elat,theta,H,p,D,n,smb,height_ratio)
% define z, the height profile

if H<1
    z = [0:1];
elseif H<100
    z = [0:1:H]; % m
else
    z = [0:10:H]; % m
end

% set parameters
c = 6; % for spherical grains (Behn et al in prep)
Tm = 273; % in K
Ts = Tm-25; % in K
dT = Tm-Ts; % K
rho = 917; % kg m^-3
cp = 2050; % J kg^-1 K^-1
K = 2.1; % W m^-1 K^-1
Acons = 2.4e-24; % Pa^-3 s^-1
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)
k0 = 11.4266; %mm^p/s (Azuma et al 2012)
k0 = k0./1000^p; %m^p/s
R = 8.314; % J/mol K
mu = 3e9; %Pa, shear modulus
%D = 0.05; %m, average grain size
g = 9.8; % m s^-2
M0 = 0.023; % m^2 s kg^-1
T = [248:1:273];
T = T(end:-1:1);

if H<100
    depth_slice = find(z==round(height_ratio*H,0));
else
    depth_slice = find(z==round(height_ratio*H,-1));
end

smb = smb./(1e3*3.154e7);
Pe = (rho.*cp.*smb.*H)./(K);

% solve for Brinkmann number
Br = ((2*theta*H^2)./(K.*(Tm-Ts)))*(elat.^(n+1)./Acons).^(1/n);

% find critical strain rate
elatcrit = ((0.5.*Pe.^2)./(Pe-1+exp(-Pe))).^(n./(n+1)).*((K.*dT)./(Acons^(-1/n).*H.^2.*theta)).^(n/(n+1));

% find the thickness of the temperate zone
if elat > elatcrit
    zetaH = 1-(Pe/(Br))-(1/Pe).*(1+real(lambertw(-exp((-Pe.^2)./(Br)-1))));
else
    zetaH = 0;
end

zetaaddiff = zetaH.*H;

% find the temperature profile
Taddiff = zeros(size(z));
for i=1:length(z)
    if and(z(i)>=zetaH.*H,z(i)<=H)
        Taddiff(i) = Ts + dT*(Br./Pe).*(1-(z(i)./H)+(1/Pe).*exp(Pe.*(zetaH-1))-(1/Pe).*exp(Pe.*(zetaH.*H-z(i))./H));
    else
        Taddiff(i) = Tm;
    end
end
temp = Taddiff(depth_slice);

% compute grain size
[grainsize_temp,tau,A] = findGrainSize(Taddiff,theta,elat,R,mu,D,M0,k0,gamma,p,c,n);
gsize = grainsize_temp(depth_slice);

% compute flow rate parameter from temperature
[Qg,Qc,Qm] = defineActivationEnergies(Taddiff);
A0 = 2.4e-25./(exp(-100./(R+273)));
A = zeros(size(z));
for i=1:length(Taddiff)
    A(i) = A0*exp(-(Qc(i)./R).*((1/Taddiff(i))-(1/263)));
end
A = A(depth_slice);

end

