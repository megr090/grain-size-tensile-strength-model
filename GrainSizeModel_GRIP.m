%% Grain Size Estimate - Grain Growth, Polygonization, Migration - Advection Diffusion

% This model estimates the grain size in a 1D model of ice frozen to its bed. 

% The steady state grain size estimate follows Austin and Evans 2007 (including a balance between 
% recrystallized grain reduction and normal grain growth) 

%% Figure 1: GRIP data, wattmeter vs my new model
clear all;

theta = 0.99;

% temperature depth profile from GRIP
load('griptemperature.mat')
tempdata = grip_temp;
tempdata(:,1) = tempdata(end:-1:1,1);

% define geometry and creep parameters
H = tempdata(1,1); % m
z = tempdata(:,1); % m
T = tempdata(:,2)+273; % Krho = 917; % kg m^-3
n = 3;
%alpha = 0.11; % degrees
alpha1 = 0.01;
alpha2 = 0.05; 
alpha1 = alpha1*(pi/180); % radians
alpha2 = alpha2*(pi/180); % radians
g = 9.8; % m s^-2
R = 8.314; % J/mol K
rho = 917; % kg m^-3

% compute viscosity
tau1 = rho.*g.*(H-z).*sin(alpha1);
tau2 = rho.*g.*(H-z).*sin(alpha2);

[Qg,Qc,Qm] = defineActivationEnergies(T); % Tc set to -10

A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));
for i=1:length(T)
    A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
end

% compute strain rate from constitutive relation
strainrate1 = A'.*tau1.^n; % 1/s
strainrate2 = A'.*tau2.^n; % 1/s

% define grain size parameters
p = 9;
c = 6; % for spherical grains (Behn et al in prep)
Tm = 273; % in K
Ts = Tm-25; % in K
dT = Tm-Ts; % K
cp = 2050; % J kg^-1 K^-1
K = 2.1; % W m^-1 K^-1
Acons = 2.4e-24; % Pa^-3 s^-1
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)
k0 = 11.4266; % mm^p/s (Azuma et al 2012)
k0 = k0./1000^p; %m^p/s
mu = 3e9; % Pa, shear modulus
D = 0.05; % m, average grain size
M0 = 0.023; % m^2 s kg^-1

% Find steady state grain size and plot depth profiles
grainsize1 = zeros(1,length(T));
grainsize2 = zeros(1,length(T));
grainsize_wattmeter = ((k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma)./(tau2.*(1-theta).*strainrate2)).^(1/(1+p));
grainsize1 = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau1.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T)))./(4.*mu.^2.*tau1.*(1-theta).*strainrate1)).^(1/(1+p));
grainsize2 = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*c.*gamma+tau2.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T)))./(4.*mu.^2.*tau2.*(1-theta).*strainrate2)).^(1/(1+p));
grainsize1 =  grainsize1.*1e3;
grainsize2 =  grainsize2.*1e3;
grainsize_wattmeter = grainsize_wattmeter.*1e3;

viscosity = zeros(1,length(T))+0.5.*(A.*strainrate1.^(n-1)).^(-1/n); % Pa s

% equilibrium time
time_to_eq1 = (p.*(1/0.05).*strainrate1).^(-1); % in s
time_to_eq1 = time_to_eq1./(365*24*60*60);
time_to_eq2 = (p.*(1/0.05).*strainrate2).^(-1); % in s
time_to_eq2 = time_to_eq2./(365*24*60*60);

premeltindx = find(tempdata(:,2) >= -13);
premeltheight = z(premeltindx(1));

% compare to GRIP data
load('gripdata.mat')
%grip_depth = grip_depth-grip_depth(end);
grip_depth = grip_depth(end:-1:1);

load('gripage.mat')
grip_depth2 = gripage(:,1);
grip_age = gripage(:,2);
%grip_depth2 = grip_depth2-grip_depth2(end);
grip_depth2 = grip_depth2(end:-1:1);

% find the depth the grain sizes are steady state
[xss1,yss1] = polyxpoly(time_to_eq1,z,grip_age,grip_depth2);
[xss2,yss2] = polyxpoly(time_to_eq2,z,grip_age,grip_depth2);

% plot profiles
figure;
subplot(1,3,1)
plot(T,z,'Color','k','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([240,270])
ylim([0 H])
hold on
plot([240:1:270],premeltheight*ones(1,length([240:1:270])),'--k','LineWidth',2)
title('Observed Temperature (K)')
%xlabel('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,3,2)
semilogx(strainrate1,z,'--b','LineWidth',2)
hold on
semilogx(strainrate2,z,'Color','b','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([10^-20 10^-10])
ylim([0 H])
plot([10^-20 10^-10],premeltheight*ones(1,length([10^-20 10^-10])),'--k','LineWidth',2)
title('Computed Strain Rate (s^{-1})')
%ylabel('Height (m)')
grid on
legend({'$$ \alpha = 0.01$$','$$ \alpha = 0.05$$'},'Interpreter','latex','location','southwest')

subplot(1,3,3)
rectangle('Position',[0 yss1 40 grip_depth(1)],'FaceColor',[0.9 .9 .9],'EdgeColor',[0.9 .9 .9])
hold on
rectangle('Position',[0 yss2 40 grip_depth(1)],'FaceColor',[0.6 .6 .6],'EdgeColor',[0.6 .6 .6])
plot(grainsize2,z,'b','LineWidth',2)
plot(grainsize_wattmeter,z,'Color','r','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
hold on
scatter(grip_size,grip_depth,'k','filled')
plot([0:1:40],premeltheight*ones(1,length([0:1:40])),'--k','LineWidth',2)
plot(grainsize1,z,'--b','LineWidth',2)
plot(grainsize2,z,'Color','b','LineWidth',2)
xlim([0,40])
ylim([0 H])
plot([0:1:40],premeltheight*ones(1,length([0:1:40])),'--k','LineWidth',2)
title('Grain Size (mm)')
%xlabel('Grain Size (mm)')
%ylabel('Height (m)') 
grid on
legend('Model, this study','Model, A+E(2007)','GRIP','location','northeast')
 