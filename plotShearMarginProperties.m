function [] = plotShearMarginProperties(SRmat,Tmat,p,icestream,Vmat,Xfix,Yfix)
Vmat(Vmat<0) = 0;

xaxis = linspace(Xfix(1,1)./1000, Xfix(1,end)./1000, size(Vmat,2));
yaxis = linspace(Yfix(1,1)./1000,Yfix(end,1)./1000,size(Vmat,1));
%Boundary = [Xfix(1,1),Xfix(end,1),Xfix(end,end),Xfix(1,end);Yfix(1,1),Yfix(end,1),Yfix(end,end),Yfix(1,end)]'./1000;

[Xa,Ya] = meshgrid(xaxis,yaxis);

Xfix = rot90(Xfix');
Yfix = rot90(Yfix');
%Boundary = [Xfix(1,1),Xfix(end,1),Xfix(end,end),Xfix(1,end);Yfix(1,1),Yfix(end,1),Yfix(end,end),Yfix(1,end)]'./1000;
Boundary = [Xfix(end,1),Xfix(1,1),Xfix(1,end),Xfix(end,end);Yfix(end,1),Yfix(1,1),Yfix(1,end),Yfix(end,end)]'./1000;

figure;
imagesc(xaxis,yaxis,real(log10(Vmat)))
colormap(colorcet('l17','reverse',0))
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xticks([]);
xticklabels({});
yticks([]);
yticklabels({});
cbar = colorbar; 
caxis([1 3.6])
cbar.Ticks=[1 2 3 3.6];
cbar.TickLabels={'10', '100', '1000', '4000'};
title('Velocity (m a^{-1})')
load('groundingline.mat')
poly = polyshape(Boundary);
for i=1:676
    lineseg2 = [];
    lineseg = [];
    lineseg = [S(i).X/1000;S(i).Y/1000];
    lineseg2(1,:) = lineseg(1,~isnan(lineseg(1,:)));
    lineseg2(2,:) = lineseg(2,~isnan(lineseg(2,:)));
    [in,out] = intersect(poly,lineseg2');
    plot(in(:,1),in(:,2),'Color','k','LineWidth',2);
    hold on;
end
set(gca,'xdir','reverse')
camroll(90)
camroll(90)
saveas(gcf,sprintf('%s_velocity.png',icestream), 'png')

figure;
poly = polyshape(Boundary);
for i=1:676
    plot(S(i).X/1000,S(i).Y/1000,'Color',[0.5 0 0],'LineWidth',3)
    hold on;
end
scatter(poly.Vertices(:,1),poly.Vertices(:,2))

SRmat = SRmat.*3.154e7;

figure;
imagesc(xaxis,yaxis,SRmat)
colormap(colorcet('l17','reverse',0))
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xticks([]);
xticklabels({});
yticks([]);
yticklabels({});
colorbar; 
caxis([0 0.15])
title('Strain Rates (a^{-1})')
load('groundingline.mat')
poly = polyshape(Boundary);
for i=1:676
    lineseg2 = [];
    lineseg = [];
    lineseg = [S(i).X/1000;S(i).Y/1000];
    lineseg2(1,:) = lineseg(1,~isnan(lineseg(1,:)));
    lineseg2(2,:) = lineseg(2,~isnan(lineseg(2,:)));
    [in,out] = intersect(poly,lineseg2');
    plot(in(:,1),in(:,2),'Color','k','LineWidth',3);
    hold on;
end
set(gca,'xdir','reverse')
camroll(90)
camroll(90)
saveas(gcf,sprintf('%s_strainrates.png',icestream), 'png')

figure;
imagesc(xaxis,yaxis,Tmat)
colormap(colorcet('l17','reverse',0))
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xticks([]);
xticklabels({});
yticks([]);
yticklabels({});
colorbar; 
caxis([0 2000])
title('Thickness (m)')
load('groundingline.mat')
poly = polyshape(Boundary);
for i=1:676
    lineseg2 = [];
    lineseg = [];
    lineseg = [S(i).X/1000;S(i).Y/1000];
    lineseg2(1,:) = lineseg(1,~isnan(lineseg(1,:)));
    lineseg2(2,:) = lineseg(2,~isnan(lineseg(2,:)));
    [in,out] = intersect(poly,lineseg2');
    plot(in(:,1),in(:,2),'Color','k','LineWidth',3);
    hold on;
end
set(gca,'xdir','reverse')
camroll(90)
camroll(90)
saveas(gcf,sprintf('%s_thickness.png',icestream), 'png')

end

