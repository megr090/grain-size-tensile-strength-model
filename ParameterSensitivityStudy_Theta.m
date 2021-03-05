%% Test the sensitivity of the grain-size model to the energy partitioning parameter Theta

%% Idealized Shear Margin Case

clear all;
H = 1000; % m
z = [0:10:H]; % m

% varying Theta
p = 9;
D = 0.05;
thetas = [0.25,0.5,0.75,0.99];
T = zeros(length(z),length(thetas));
grainsize = zeros(length(z),length(thetas));
elat = 1e-9;
n = 3;
for theta = 1:length(thetas)
    [T(:,theta),grainsize(:,theta),zetaaddiff] = computeGrainSizefromSR(elat,D,p,n,thetas(theta),H,z);
end

cmap = colorcet('l16');
cmap = cmap(1:52:end,:);
cmap = cmap(2:end,:);
cmap = cmap(end:-1:1,:);

% plot temperature profiles
figure;
subplot(1,2,1)
for i=1:length(thetas)
    plot(T(:,i),z,'Color',cmap(i,:),'LineWidth',2)
    hold on
end
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([273-25,273])
ylim([0 H])
title('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,2,2)
for i=1:length(thetas)
    plot(grainsize(:,i),z,'LineWidth',2,'Color',cmap(i,:))
    hold on
end
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
title('Grain Size (mm)')
grid on
legend({'$$ \Theta = 0.25 $$','$$ \Theta = 0.5 $$','$$ \Theta = 0.75 $$','$$ \Theta \approx 1 $$'},'Interpreter','latex')

%% Function
function [Taddiff,grainsize_temp,zetaaddiff] = computeGrainSizefromSR(elat,D,p,n,theta,H,z)
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
f = theta;
Pe = 2;

% solve for Brinkmann number
Br = ((2*f*H^2)./(K.*(Tm-Ts)))*(elat.^(n+1)./Acons).^(1/n);

% find critical strain rate
elatcrit = ((0.5.*Pe.^2)./(Pe-1+exp(-Pe))).^(n./(n+1)).*((K.*dT)./(Acons^(-1/n).*H.^2.*f)).^(n/(n+1));

% find the thickness of the temperate zone
if elat > elatcrit
    zetaH = 1-(Pe/(Br))-(1/Pe).*(1+real(lambertw(-exp((-Pe.^2)./(Br)-1))));
else
    zetaH = 0;
end

zetaaddiff = zetaH.*H;

% find the temperature profile
for i=1:length(z)
    if and(z(i)>=zetaH.*H,z(i)<=H)
        Taddiff(i) = Ts + dT*(Br./Pe).*(1-(z(i)./H)+(1/Pe).*exp(Pe.*(zetaH-1))-(1/Pe).*exp(Pe.*(zetaH.*H-z(i))./H));
    else
        Taddiff(i) = Tm;
    end
end

[grainsize_temp,tau,A] = findGrainSize(Taddiff,f,elat,R,mu,D,M0,k0,gamma,p,c,n);

end