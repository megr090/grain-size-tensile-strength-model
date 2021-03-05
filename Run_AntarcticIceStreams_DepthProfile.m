%% Create plots from Theta of Antarctic Ice Streams

clear all;
% choose the ice stream
icestream = 'PIG'; % options: BM (Bindschadler/MacAyeal), Byrd, Recovery, Amery, PIG

[SRmat,Tmat,Vmat,SMBmat,Xfix,Yfix] = readData(icestream);

% define key parameters
n = 3;
p = 2;
D = 0.05;
% Br = 4;
% Pe = 2;
% elat = 1e-9;

festmat = loadThetaEstimate(icestream,p);
festmat = 0.99.*ones(size(festmat));

% Find Shear Margin Properties: Viscosity, A, Temperate Ice Thickness,
% Grain Reduction, Grain Size for various depths
i = 324;
j = 119;
[viscositymat,Amat,zetamat,redmat,gsizemat,tmpmat,z,H] = findShearMarginProperties_Antarctica_DProf(SRmat,festmat,Tmat,SMBmat,p,D,n,i,j);

% Tensile Strength Depth Slices
K = 0.052; % MPa m^(1/2)
tensstrength_depthprof = K.*(gsizemat./(1e3)).^(-0.5);

figure;
plot(tmpmat.*(14./25)-((273-25)*14/25),z,'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880])
set(gca,'FontSize',14,'FontWeight','b','GridColor','r');
grid on
ylim([0 H])
ylabel('Height (m)')
yticks([H/4,H/2,(3*H)/4,H])
yticklabels({'0.25H','0.5H','0.75H','H=427'})
ax1 = gca; % current axes
ax1.XLim = [0 14];
ax1.XColor = [0.8500, 0.3250, 0.0980];
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax1.Position = ax1_pos + [0 0.3 0 -0.3];
hold on;
ax1.XLabel.String = 'Grain Size (mm)';
%ax2 = axes('Position',ax1_pos,...
%    'XAxisLocation','top',...
%    'YAxisLocation','right',...
%    'Color','none');
ax2 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.15 0 0], 'color', 'none', 'linewidth', 2);
ax2.XColor = [0, 0.4470, 0.7410];
ax2.YColor = 'k';
line(tensstrength_depthprof.*14./8,z,'LineWidth',2,'Parent',ax1,'Color',[0, 0.4470, 0.7410])
ax2.XLim = [0 8];
ax2.XLabel.String = 'Tensile Strength (MPa)';
set(ax2,'FontSize',14,'FontWeight','b','GridColor','r');
ax2.YAxis.Visible = 'off';
ax3 = axes('position', (ax1_pos .* [1 1 1 1e-3]) + [0 0.00 0 0], 'color', 'none', 'linewidth', 2);
line(gsizemat,z,'LineWidth',2,'Parent',ax1,'Color',[0.8500, 0.3250, 0.0980])
ax3.XLim = [273-25 273];
ax3.XLabel.String = 'Ice Temperature (K)';
ax3.XColor = [0.4660, 0.6740, 0.1880];
ax3.YColor = 'b';
set(ax3,'FontSize',14,'FontWeight','b','GridColor','r');
ax3.YAxis.Visible = 'off';
