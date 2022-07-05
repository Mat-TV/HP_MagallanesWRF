%Función para graficar un perfil vertical latitudinal  
  % dentro de la función están las variables PATH{in,out}
%Argumentos de entrada:
    %filename: nombre del archivo
%Argumentos de salida:
    %ninguno: nada
%  2022-06-25 - Mat Troncoso Villar (matro1432@gmail.com)
function perfil_theta(filename)
% aseo
% filename='wrfout_d03_2020-01-03';
PATHin = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\exp006\';
PATHout = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\Plots\';
lat=ncread([PATHin,filename],'XLAT');
lat = double(lat(:,:,1));
%
lon=ncread([PATHin,filename],'XLONG');
lon = double(lon(:,:,1));
%
lon_plot = lon(:,1);
lat_plot = lat(1,:)';
%
xx = -67.61; % longitud a la que realizar el perfil
[minimox,positionx]=min(abs(lon_plot-xx));
mas_cercanox=lon(positionx);
%
hgt=ncread([PATHin,filename],'HGT',[1,1,1],[100,100,1]);
hgt = double(hgt); % (longitud,latitud)
hgt = hgt(positionx,:);
wwind=ncread([PATHin,filename],'W',[1,1,1,1],[100,100,44,121]);
wwind = squeeze(double(wwind(positionx,:,:,:))); % (latitud,verticales,t)
wwind = mean(wwind,3).*100; % velocidad vertical [cm/s]
wwind = wwind(:,2:end);
%
theta = ncread([PATHin,filename],'T',[1,1,1,1],[100,100,43,121]);%
theta = squeeze(double(theta(positionx,:,:,:))); % (latitud,verticales,t)
theta = theta + 300 -273.15; % Temperatura potencial [°C]
theta = mean(theta,3); % 
%
P=ncread([PATHin,filename],'P',[1 1 1 1],[1 1 Inf 1]);
P=squeeze(double(P));
PB=ncread([PATHin,filename],'PB',[1 1 1 1],[1 1 Inf 1]);
PB=squeeze(double(PB));
PT=(P+PB)./100;% presion total [hPa]
clear minimox P PB 
%
%% convertir la altura a presión
% Considerando atm isotérmica T=cte
R=287; %constante de gas de aire seco [J/kgK]
T0=237.1; %temperatura promedio [K]
g=9.8; % [m/s^2]
%
H= R*T0/g; %escala de altura [m]
%
z=hgt; %altura [m]
p0=1013.23; %presión al nivel del mar [hPa]
p=p0*exp(-z/H);
supP = zeros(size(p))+p0;
%
%%
wwind_pos = max(wwind,5);
wwind_pos(wwind_pos==5) = NaN;
wwind_pos(wwind_pos>40) = NaN;
wwind_neg = min(wwind,-5);
wwind_neg(wwind_neg==-5) = NaN;
wwind_neg(wwind_neg<-40) = NaN;
pintada = gray(100);
%
figure
ax1=axes;
ccpos=contourf(lat_plot,PT,wwind_pos');
% colormap(pintada(51:100,:))
hold on
ccneg=contourf(lat_plot,PT,abs(wwind_neg)','--');
colormap(pintada(51:100,:))
% shading interp
ylabel('Presión [hPa]')
xlabel('Latitud [°]')
title("Sección meridional a " + mas_cercanox + " W, \theta [°C], y viento vertical [cm/s]")
ax2=axes;
cct=contour(lat_plot,PT,theta',[15 25 35 50 75 100 200],'k','ShowText','on','linewidth',1.5);
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
hf = fill_between(lat_plot,p,supP);
hf.FaceColor = [.01 .01 .01];
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.88 .11 .0675 .815]);
cb1.Label.String = 'Velocidad vertical [cm/s] (negativo -> segmentado)';
set([ax1 ax2], 'YDir','reverse')
hold off
  set(gcf, 'PaperUnits', 'centimeters')
  set(gcf, 'PaperPosition', [0 0 20 15]);
  saveas(gcf,[PATHout,'perfil_theta'],'png')
end