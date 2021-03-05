%% Grain Size Estimate - Grain Growth, Polygonization, Migration 

% This model estimates the grain size in a 1D thermomechanical model that
% considers advection-diffusion of ice.  

% The steady state grain size estimate follows Austin and Evans 2007 (including a balance between 
% recrystallized grain reduction and normal grain growth) and adds a
% parameterization of migration recrystallization

% Meghana Ranganathan [meghanar@mit.edu]

%% Large-Scale Parameters

clear all;

% define parameters for models
H = 1000; % ice thickness, m
z = [0:10:H]; % m
N = 1000./10;
n = 3; % flow law exponent
p = 2; % grain-growth exponent
c = 6; % for spherical grains (Behn et al in prep)
Tm = 273; % melting temperature, in K
Ts = Tm-25; % surface temperature, in K
dT = Tm-Ts; % K
rho = 917; % ice density, kg m^-3
cp = 2050; % J kg^-1 K^-1
K = 2.1; % W m^-1 K^-1
Acons = 2.4e-24; % prefactor in Glen's flow law, Pa^-3 s^-1
gamma = 0.065; % grain boundary energy, J/m^2 (Cuffey and Paterson 2010)
k0 = 11.4266; %mm^p/s (Azuma et al 2012)
k0 = k0./1000^p; % prefactor in the Arrhenius relation for grain growth, m^p/s
R = 8.314; % J/mol K
mu = 3e9; %Pa, shear modulus
D = 0.05; %m, average grain size
g = 9.8; % m s^-2
M0 = 0.023; % intrinsic grain boundary mobility, m^2 s kg^-1
theta = 0.99; % fraction of energy going into heating

%% Compute Ice Temperature and Grain Size

% input the Brinkmann number and Peclet number
Br = theta.*4; % scale the Brinkmann number by fraction of energy going into heating
Pe = 2; % Peclet number
 
[Taddiff,zetaaddiff,elat,elatcrit] = findIceTemperature(Br,Pe,theta,H,z,Tm,Ts,K,Acons,n,dT);

[grainsize_temp,tau,A] = findGrainSize(Taddiff,theta,elat,R,mu,D,M0,k0,gamma,p,c,n);


%% Make plots

[done] = plotGrainSizeIceTemperature(elat,z,H,Taddiff,grainsize_temp,tau,p,zetaaddiff,N);
