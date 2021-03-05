%% Compute Grain Size and Tensile Strength in Idealized Setting
% This code computes grain size and tensile strength from a set strain rate
% and plots Figure 3b in the submitted paper.

% Meghana Ranganathan [meghanar@mit.edu]

clear all;

% define f, D, p
theta = 0.99;
D = 0.05;
p = 2;

% define strain rates:
elats = [6e-10,1.3e-9,6e-9];
for elat = 1:length(elats)
    [Taddiff,grainsize_temp,zetaaddiff,z,H] = findtempgsfromenergy(elats(elat),D,p,theta);
    grainsize(:,elat) = grainsize_temp;
    zetaaddiff_t(elat) = zetaaddiff;
    Taddiff_t(:,elat) = Taddiff;
end

zetaaddiff_t(zetaaddiff_t==0) = NaN;

cmap = colorcet('fire');
cmap = cmap(1:86:end,:);

figure;
hold on;
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([273-25,273])
ylim([0 H])
plot(Taddiff_t(:,3),zetaaddiff_t(3).*ones(length(z),1),'--k','LineWidth',2)
inBetween = [z, fliplr(z)];
x2 = [Taddiff_t(:,1)', fliplr(Taddiff_t(:,3)')];
plot(Taddiff_t(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(Taddiff_t(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(Taddiff_t(:,1),z,'LineWidth',3,'Color',cmap(1,:));
fill(x2,inBetween, [0.9 0.9 0.9]);
plot(Taddiff_t(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(Taddiff_t(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(Taddiff_t(:,1),z,'LineWidth',3,'Color',cmap(1,:));
plot(Taddiff_t(:,3),zetaaddiff_t(3).*ones(length(z),1),'--','Color',cmap(3,:),'LineWidth',2)
plot(Taddiff_t(:,2),zetaaddiff_t(2).*ones(length(z),1),'--','Color',cmap(2,:),'LineWidth',2)
plot(Taddiff_t(:,1),zetaaddiff_t(1).*ones(length(z),1),'--','Color',cmap(1,:),'LineWidth',2)
ylabel('Height Above Bed (m)')
xlabel('Temperature (K)')
legend('Temperate Zone')

K = 0.052; % MPa m^(1/2)
kt = 0.03; % MPa m^(1/2)
sigma0 = 0.52; % MPa
dc = ((K-kt)/sigma0).^2;

grainsize_c = grainsize;
grainsize = log(grainsize);

figure;
hold on;
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
plot(grainsize(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(grainsize(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(grainsize(:,1),z,'LineWidth',3,'Color',cmap(1,:));
inBetween = [z, fliplr(z)];
x2 = [grainsize(:,1)', fliplr(grainsize(:,3)')];
fill(x2,inBetween, [0.9 0.9 0.9]);
plot(grainsize(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(grainsize(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(grainsize(:,1),z,'LineWidth',3,'Color',cmap(1,:));
xticks([-3 -2 -1 0 1 2 3])
xlabel('Grain Size (mm)')
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2','10^3'})
%xlim([0 60])
legend({'$$ \dot{\epsilon} = 6 \times 10^{-9} s^{-1}$$','$$ \dot{\epsilon} = 1.3 \times 10^{-9} s^{-1}$$','$$ \dot{\epsilon} = 6 \times 10^{-10} s^{-1}$$'},'Interpreter','latex','location','southeast');

gsize_m = grainsize_c./(1e3);
K = 0.052; % MPa m^(1/2)
tens_strength_prop = K.*gsize_m.^(-0.5);
tens_strength_prop = log(tens_strength_prop);

figure;
hold on;
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylim([0 H])
plot(tens_strength_prop(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(tens_strength_prop(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(tens_strength_prop(:,1),z,'LineWidth',3,'Color',cmap(1,:));
inBetween = [z, fliplr(z)];
x2 = [tens_strength_prop(:,1)', fliplr(tens_strength_prop(:,3)')];
fill(x2,inBetween, [0.9 0.9 0.9]);
plot(tens_strength_prop(:,3),z,'LineWidth',3,'Color',cmap(3,:));
plot(tens_strength_prop(:,2),z,'LineWidth',3,'Color',cmap(2,:));
plot(tens_strength_prop(:,1),z,'LineWidth',3,'Color',cmap(1,:));
xlabel('Tensile Strength (MPa)')
xticks([-2 -1 0 1 2])
xticklabels({'10^{-2}','10^{-1}','10^0','10^1','10^2'})
%xlim([0 45])

%%
function [Taddiff,grainsize_temp,zetaaddiff,z,H] = findtempgsfromenergy(elat,D,p,theta)
H = 1000; % m
z = [0:10:H]; % m
N = H./10;
n = 3;
%p = 9;
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

Pe = 2;

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
for i=1:length(z)
    if and(z(i)>=zetaH.*H,z(i)<=H)
        Taddiff(i) = Ts + dT*(Br./Pe).*(1-(z(i)./H)+(1/Pe).*exp(Pe.*(zetaH-1))-(1/Pe).*exp(Pe.*(zetaH.*H-z(i))./H));
    else
        Taddiff(i) = Tm;
    end
end

[grainsize_temp,tau,A] = findGrainSize(Taddiff,theta,elat,R,mu,D,M0,k0,gamma,p,c,n);

end
