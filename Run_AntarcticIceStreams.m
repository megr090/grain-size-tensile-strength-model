%% Compute Grain Size, Ice Temperature in Antarctic Ice Streams
% Here we compute grain size, ice temperature, and other properties of
% shear margins in 5 Antarctic ice streams: Bindschadler/MacAyeal Ice
% Streams, Byrd Glacier, Recovery Glacier, Amery, and Pine Island Glacier.
% The necessary data is: surface mass balance, strain rates, ice thickness,
% surface velocity

% Meghana Ranganathan [meghanar@mit.edu]

clear all;

% choose the ice stream
icestream = 'PIG'; % options: BM (Bindschadler/MacAyeal), Byrd, Recovery, Amery, PIG

% read data
[SRmat,Tmat,Vmat,SMBmat,Xfix,Yfix] = readData(icestream);

% define key parameters
n = 3;
p = 2;
D = 0.05;
theta = 0.99;
height_ratio = 0.75; % 0 0.25 0.5 0.75

% Find Shear Margin Properties: Viscosity, A, Temperate Ice Thickness,
% Grain Reduction, Grain Size for various depths
[Amat,gsizemat,tmpmat] = findShearMarginProperties_Antarctica(SRmat,theta,Tmat,SMBmat,p,D,n,height_ratio);
