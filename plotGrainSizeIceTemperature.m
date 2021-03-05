function [done] = plotGrainSizeIceTemperature(elat,z,H,Taddiff,grainsize_temp,tau,p,zetaaddiff,N)

% compute equilibrium time
time_to_eq = (p.*(1/0.05).*elat).^(-1); % in s
time_to_eq = time_to_eq./(365*24*60*60);

% compute age depth
accumulation = (elat.*365.*24.*60.*60).*H;
age = (H./accumulation).*log(1./(1-(z./H)));
age = age(end:-1:1);

%z=z-H;

% plot time to equilibrium and age-depth
figure;
plot(age,z,'r','LineWidth',2)
hold on
%xlim([0 2.5e5])
ylim([0 H])
plot(time_to_eq+zeros(size(age)),z,'b','LineWidth',2)
xlabel('Time (yr)')
ylabel('Depth (m)')
legend('Age','Time to Equilibrium')
grid on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');


% plot temperature profiles
figure;
subplot(1,3,1)
%plot(Taddiff,z,'Color','r','LineWidth',2)
surface([Taddiff;Taddiff],[z;z],[zeros(size(z));zeros(size(z))],[Taddiff;Taddiff],...
        'FaceColor','none',...
        'EdgeColor','interp','LineWidth',3);
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([273-25,273])
ylim([0 H])
plot(Taddiff,zetaaddiff*ones(1,N+1),'--r','LineWidth',2)
title('Temperature (K)')
%xlabel('Temperature (K)')
ylabel('Height (m)')
grid on

subplot(1,3,2)
%plot(tau.*1e-6,z,'Color','r','LineWidth',2)
surface([tau.*1e-6;tau.*1e-6],[z;z],[zeros(size(z));zeros(size(z))],[Taddiff;Taddiff],...
        'FaceColor','none',...
        'EdgeColor','interp','LineWidth',3);
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
xlim([0 0.5])
ylim([0 H])
title('Shear Stress (MPa)')
%xlabel('Shear Stress (MPa)')
%ylabel('Height (m)')
grid on

subplot(1,3,3)
%plot(grainsize_temp,z,'Color','r','LineWidth',2)
surface([grainsize_temp;grainsize_temp],[z;z],[zeros(size(z));zeros(size(z))],[Taddiff;Taddiff],...
        'FaceColor','none',...
        'EdgeColor','interp','LineWidth',3);
hold on
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
%plot(grainsize_temp,zetaaddiff*ones(1,N+1),'--r')
%xlim([273-25,273])
ylim([0 H])
title('Grain Size (mm)')
%xlabel('Grain Size (mm)')
%ylabel('Height (m)') 
grid on

done = 1;
end

