%% Test the sensitivity of the grain-size model to the critical temperature

%% Idealized Shear Margin Case

clear all;
H = 1000; % m
z = [0:10:H]; % m

% varying Q's
p = 9;
D = 0.05;
theta = 0.99;
elat = 1e-9;
n = 3;
tcs = [-13:1:-8];
Qgm = 40;
Qgp = 100;

grainsize = zeros(length(z),length(tcs));
for tc = 1:length(tcs)
    [T,grainsize(:,tc),zetaaddiff] = computeGrainSizefromSR(elat,D,p,n,theta,H,z,Qgm,Qgp,tcs(tc));
end

cmap = colorcet('l16');
cmap = cmap(1:36:end,:);
cmap = cmap(2:end,:);
cmap = cmap(end:-1:1,:);

% plot temperature profiles
figure;
subplot(1,2,1)
plot(T,z,'Color','k','LineWidth',2)
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([273-25,273])
ylim([0 H])
title('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,2,2)
for i=1:length(tcs)
    plot(grainsize(:,i),z,'LineWidth',2,'Color',cmap(i,:))
    hold on
end
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
title('Grain Size (mm)')
grid on
legend({'$$ T_c = -13 $$','$$ T_c = -12 $$','$$ T_c = -11 $$','$$ T_c = -10 $$','$$ T_c = -9 $$','$$ T_c = -8 $$'},'Interpreter','latex')

%% Comparisons to GRIP Data

clear all;

% temperature depth profile from GRIP
load('griptemperature.mat')
tempdata = grip_temp;
tempdata(:,1) = tempdata(end:-1:1,1);

% define geometry and creep parameters
H = tempdata(1,1); % m
z = tempdata(:,1); % m
T = tempdata(:,2)+273; % K
rho = 917; % kg m^-3
alpha = 0.05; % degrees
alpha = alpha*(pi/180); % radians
g = 9.8; % m s^-2
R = 8.314; % J/mol K
A0 = 2.4e-25./(exp(-100./(R+273)));
n = 3;
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)

% varying tc's
tcs = [-13:1:-8];
Qgm = 40;
Qgp = 100;

% compute stress
tau = rho.*g.*(H-z).*sin(alpha);

c = 6;
f = 0.99;
mu = 3e9;
M0 = 0.023; % m^2 s kg^-1
p = 9;
k0 = 11.4266; %mm^p/s (Azuma et al 2012) at Q=42 kJ/mol
k0 = k0./1000.^p; %m^p/s
D = 0.05;
for tc=1:length(tcs)
    [Qg,Qc,Qm] = defineActivationEnergies_Test(T,Qgm,Qgp,tcs(tc));
    
    % compute A
    A = zeros(size(T));
    for i=1:length(T)
        A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
    end
    
    % compute strain rate from constitutive relation
    strainrate = A.*tau.^n; % 1/s
    
    [grainsize(:,tc),viscosity,time_to_eq,tau] = findGrainSize_TestAE(T,f,strainrate,R,mu,D,M0,k0,gamma,p,c,n,Qgm,Qgp,tcs(tc));
end

% equilibrium time
time_to_eq = (p.*(1/0.05).*strainrate).^(-1); % in s
time_to_eq = time_to_eq./(365*24*60*60);

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
[xss,yss] = polyxpoly(time_to_eq,z,grip_age,grip_depth2);

cmap = colorcet('l16');
cmap = cmap(1:36:end,:);
cmap = cmap(2:end,:);
cmap = cmap(end:-1:1,:);

% plot temperature profiles
figure;
subplot(1,2,1)
plot(T,z,'Color','r','LineWidth',2)
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([240,270])
ylim([0 H])
title('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,2,2)
rectangle('Position',[0 yss 40 grip_depth(1)],'FaceColor',[0.8 .8 .8],'EdgeColor',[0.8 .8 .8])
hold on
for i=1:length(tcs)
    plot(grainsize(:,i),z,'LineWidth',2,'Color',cmap(i,:))
    hold on
end
scatter(grip_size,grip_depth,'k')
xlim([0 40])
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
title('Grain Size (mm)')
grid on
legend({'$$ T_c = -13 $$','$$ T_c = -12 $$','$$ T_c = -11 $$','$$ T_c = -10 $$','$$ T_c = -9 $$','$$ T_c = -8 $$','GRIP'},'Interpreter','latex')

%% Comparisons to GRIP Data

clear all;

% temperature depth profile from GISP2
load('gisp2temperature.mat')
tempdata = gisp2_temp;
tempdata(:,1) = tempdata(end:-1:1,1);

% define geometry and creep parameters
H = tempdata(1,1); % m
z = tempdata(:,1); % m
T = tempdata(:,2)+273; % K
rho = 917; % kg m^-3
alpha = 0.05; % degrees
alpha = alpha*(pi/180); % radians
g = 9.8; % m s^-2
R = 8.314; % J/mol K
A0 = 2.4e-25./(exp(-100./(R+273)));
n = 3;
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)

% varying tc's
tcs = [-13:1:-8];
Qgm = 40;
Qgp = 100;

% compute stress
tau = rho.*g.*(H-z).*sin(alpha);

c = 6;
f = 0.99;
mu = 3e9;
M0 = 0.023; % m^2 s kg^-1
p = 9;
k0 = 11.4266; %mm^p/s (Azuma et al 2012) at Q=42 kJ/mol
k0 = k0./1000.^p; %m^p/s
D = 0.05;
for tc=1:length(tcs)
    [Qg,Qc,Qm] = defineActivationEnergies_Test(T,Qgm,Qgp,tcs(tc));
    
    % compute A
    A = zeros(size(T));
    for i=1:length(T)
        A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
    end
    
    % compute strain rate from constitutive relation
    strainrate = A.*tau.^n; % 1/s
    
    [grainsize_model(:,tc),viscosity,time_to_eq,tau] = findGrainSize_TestAE(T,f,strainrate,R,mu,D,M0,k0,gamma,p,c,n,Qgm,Qgp,tcs(tc));
end

% equilibrium time
time_to_eq = (p.*(1/0.05).*strainrate).^(-1); % in s
time_to_eq = time_to_eq./(365*24*60*60);

% compare to GISP2 data
load('gisp2age.mat')
gisp2_depth2 = gisp2age(:,1);
gisp2_age = gisp2age(:,2);
%gisp2_depth2 = gisp2_depth2-gisp2_depth2(end);
gisp2_depth2 = gisp2_depth2(end:-1:1);

%z = z-H;

load('gisp2grainsize.mat')
gisp2_depth = depth(end:-1:1);
gisp2_depth = gisp2_depth.*1e3;
gisp2_grainsize = grainsize.*1e3;

% find the depth the grain sizes are steady state
[xss,yss] = polyxpoly(time_to_eq,z,gisp2_age,gisp2_depth2);

cmap = colorcet('l16');
cmap = cmap(1:36:end,:);
cmap = cmap(2:end,:);
cmap = cmap(end:-1:1,:);

% plot temperature profiles
figure;
subplot(1,2,1)
plot(T,z,'Color','r','LineWidth',2)
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([240,270])
ylim([0 H])
title('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,2,2)
rectangle('Position',[0 yss 40 gisp2_depth(1)],'FaceColor',[0.8 .8 .8],'EdgeColor',[0.8 .8 .8])
hold on
for i=1:length(tcs)
    plot(grainsize_model(:,i),z,'LineWidth',2,'Color',cmap(i,:))
    hold on
end
scatter(gisp2_grainsize,gisp2_depth,'k')
xlim([0 40])
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
title('Grain Size (mm)')
grid on
legend({'$$ T_c = -13 $$','$$ T_c = -12 $$','$$ T_c = -11 $$','$$ T_c = -10 $$','$$ T_c = -9 $$','$$ T_c = -8 $$','GISP2'},'Interpreter','latex')

%% Comparisons to WAIS Divide Data

clear all;

% temperature depth profile from GRIP
load('waisdividetemperature.mat')
tempdata = waisdivide_temp(1:100:end,:);
tempdata(:,1) = tempdata(end:-1:1,1);

% define geometry and creep parameters
H = tempdata(1,1); % m
z = tempdata(:,1); % m
T = tempdata(:,2)+273; % K
rho = 917; % kg m^-3
alpha = 0.05; % degrees
alpha = alpha*(pi/180); % radians
g = 9.8; % m s^-2
R = 8.314; % J/mol K
A0 = 2.4e-25./(exp(-100./(R+273)));
n = 3;
gamma = 0.065; % J/m^2 (Cuffey and Paterson 2010)

% varying tc's
tcs = [-13:1:-8];
Qgm = 40;
Qgp = 100;
grainsize = zeros(length(z),length(tcs));

% compute viscosity
tau = rho.*g.*(H-z).*sin(alpha);

c = 6;
f = 0.99;
mu = 3e9;
M0 = 0.023; % m^2 s kg^-1
p = 9;
k0 = 11.4266; %mm^p/s (Azuma et al 2012) at Q=42 kJ/mol
k0 = k0./1000.^p; %m^p/s
D = 0.05;
for tc=1:length(tcs)
    [Qg,Qc,Qm] = defineActivationEnergies_Test(T,Qgm,Qgp,tcs(tc));
    
    % compute A
    A = zeros(size(T));
    for i=1:length(T)
        A(i) = A0*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
    end
    
    % compute strain rate from constitutive relation
    strainrate = A.*tau.^n; % 1/s
    
    [grainsize(:,tc),viscosity,time_to_eq,tau] = findGrainSize_TestAE(T,f,strainrate,R,mu,D,M0,k0,gamma,p,c,n,Qgm,Qgp,tcs(tc));
end

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

cmap = colorcet('l16');
cmap = cmap(1:36:end,:);
cmap = cmap(2:end,:);
cmap = cmap(end:-1:1,:);

% plot temperature profiles
figure;
subplot(1,2,1)
plot(T,z,'Color','r','LineWidth',2)
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([240,270])
ylim([0 H])
title('Temperature (K)')
ylabel('Height Above Bed (m)')
grid on

subplot(1,2,2)
rectangle('Position',[0 yss 40 waisdivide_depth(1)],'FaceColor',[0.8 .8 .8],'EdgeColor',[0.8 .8 .8])
hold on
for i=1:length(tcs)
    plot(grainsize(:,i),z,'LineWidth',2,'Color',cmap(i,:))
    hold on
end
scatter(waisdivide_size,waisdivide_depth,'k')
xlim([0 40])
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
title('Grain Size (mm)')
grid on
legend({'$$ T_c = -13 $$','$$ T_c = -12 $$','$$ T_c = -11 $$','$$ T_c = -10 $$','$$ T_c = -9 $$','$$ T_c = -8 $$','WAIS Divide'},'Interpreter','latex')

%% Function
function [Taddiff,grainsize_temp,zetaaddiff] = computeGrainSizefromSR(elat,D,p,n,theta,H,z,Qgm,Qgp,tc)
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

[grainsize_temp,viscosity,time_to_eq,tau] = findGrainSize_TestAE(Taddiff,f,elat,R,mu,D,M0,k0,gamma,p,c,n,Qgm,Qgp,tc);

end