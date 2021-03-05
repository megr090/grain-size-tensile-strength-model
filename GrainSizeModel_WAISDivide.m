%% Grain Size Estimate - Grain Growth, Polygonization, Migration - Advection Diffusion

% This model estimates the grain size in a 1D model of ice frozen to its bed. 

% The steady state grain size estimate follows Austin and Evans 2007 (including a balance between 
% recrystallized grain reduction and normal grain growth) 

%% Large-Scale Parameters

clear all;
theta = 0.99; % fraction of energy going into heating

%% Temperature Profile

% temperature depth profile from GRIP
load('waisdividetemperature.mat')
tempdata = waisdivide_temp(1:100:end,:);
tempdata(:,1) = tempdata(end:-1:1,1);

% define geometry and creep parameters
H = tempdata(1,1); % m
z = tempdata(:,1); % m
T = tempdata(:,2)+273; % K
rho = 917; % kg m^-3
n = 3;
alpha = 0.05; % degrees
alpha = alpha*(pi/180); % radians
g = 9.8; % m s^-2
R = 8.314; % J/mol K
A0 = 2.4e-25./(exp(-100./(R+273)));

[Qg,Qc,Qm] = defineActivationEnergies(T);

% compute viscosity
tau = rho.*g.*(H-z).*sin(alpha);

% compute A
for i=1:length(T)
    A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
end

% compute strain rate from constitutive relation
strainrate = A'.*tau.^n; % 1/s

%% Find steady state grain size and plot depth profiles

% define grain size parameters
p = 9;
cons1 = 6; % for spherical grains (Behn et al in prep)
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)
k0 = 11.4266; %mm^p/s (Azuma et al 2012)
k0 = k0./1000^p; %m^p/s
n = 3;
R = 8.314; % J/mol K
mu = 3e9; %Pa, shear modulus
D = 0.05; %m, average grain size
M0 = 0.023; % m^2 s kg^-1

grainsize_temp = zeros(length(T),1);
grainsize_temp = ((4.*mu.^2.*k0.*exp(-Qg./(R.*T)).*p.^(-1).*cons1.*gamma+tau.^4.*D.^(p).*(0.5.*p).*M0.*exp(-Qm./(R.*T)))./(4.*mu.^2.*tau.*(1-theta).*strainrate)).^(1/(1+p));
grainsize_temp =  grainsize_temp.*1e3;

% equilibrium time
time_to_eq = (p.*(1/0.05).*strainrate).^(-1); % in s
time_to_eq = time_to_eq./(365*24*60*60);

% compare to GRIP data
load('waisdividegrainsize.mat')
waisdivide_depth = waisdivide_depth(end:-1:1);

load('waisdivideage.mat')
waisdivide_depth2 = waisdivideage(:,1);
waisdivide_age = waisdivideage(:,2)*1000;
%waisdivide_depth2 = waisdivide_depth2-waisdivide_depth2(end);
waisdivide_depth2 = waisdivide_depth2(end:-1:1);

%z = z-H;

%find the depth the grain sizes are steady state
[xss,yss] = polyxpoly(time_to_eq,z,waisdivide_age,waisdivide_depth2);

% plot time to equilibrium and age-depth
figure;
plot(waisdivide_age,waisdivide_depth2,'r','LineWidth',2)
hold on
xlim([0 2.5e5])
ylim([0 H])
plot(time_to_eq,z,'b','LineWidth',2)
xlabel('Time (yr)')
ylabel('Depth (m)')
legend('WAIS Divide Age','Time to Equilibrium')
grid on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');

% plot profiles
figure;
subplot(1,4,1)
plot(T,z,'Color','k','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([240,270])
ylim([0 H])
hold on
title('Observed Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,4,2)
semilogx(tau.*1e-6,z,'b','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([10^-4 10^0])
ylim([0 H])
title('Computed Shear Stress (MPa)')
grid on

subplot(1,4,3)
semilogx(strainrate,z,'b','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([10^-20 10^-10])
ylim([0 H])
title('Computed Strain Rate (s^{-1})')
grid on

subplot(1,4,4)
rectangle('Position',[0 yss 40 waisdivide_depth(1)],'FaceColor',[0.8 .8 .8],'EdgeColor',[0.8 .8 .8])
hold on
plot(grainsize_temp,z,'Color','b','LineWidth',2)
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
hold on
scatter(waisdivide_size,waisdivide_depth,'k')
xlim([0,40])
ylim([0 H])
title('Grain Size (mm)')
%xlabel('Grain Size (mm)')
%ylabel('Height (m)') 
grid on
legend('Model, This Study','WAIS Divide','location','northeast')
