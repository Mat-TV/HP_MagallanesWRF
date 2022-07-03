%Función para graficar la topografía de un dominio de WRFOUT 
  % dentro de la función están las variables PATH{in,out}
%Argumentos de entrada:
    %filename: nombre del archivo
%Argumentos de salida:
    %ninguno: nada
%  2022-06-25 - Mat Troncoso Villar (matro1432@gmail.com)
function mapa_topografia(filename)
% aseo
% filename='wrfout_d03_2020-01-03';
PATHin = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\exp006\';
PATHout = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\Plots\';
lat=ncread([PATHin,filename],'XLAT');
lat = lat(:,:,1);
lon=ncread([PATHin,filename],'XLONG');
lon = lon(:,:,1);
hgt=ncread([PATHin,filename],'HGT',[1,1,1],[100,100,1]);
%
%%
figure
m_proj('lambert','lon',[min(min(double(lon))) max(max(double(lon)))],...
        'lat',[min(min(double(lat))) max(max(double(lat)))]);
m_pcolor(lon(:,:),lat(:,:),hgt(:,:));
m_grid('box','fancy','tickdir','in');
colormap(flipud(gray));
ax1=m_contfbar(.97,[.5 .9],hgt(:,:), ...
    [64],'edgecolor','none','endpiece','no');
title(ax1,{'Altitud [m]',''}); % Move up by inserting a blank line
% colormap([m_colmap('blues',80);m_colmap('gland',48)]);
% m_gshhs_h('color','k'); % por si es que no se guarda gumby antes
m_usercoast('gumby','linewidth',1.2,'color','k'); %alternativa a m_gshhs_h
title('Topografía del dominio 3.')
xlabel('Longitud','fontsize',12)
ylabel('Latitud','fontsize',12)
%   x_width=17.12 ;y_width=10.34;
%   cd 'C:\Users\fredo\OneDrive\Documents\Geof\Modelación numérica de la atmósfera\Trabajo 1\Gráficos';
%   set(gcf, 'PaperPosition', [0 0 x_width y_width]);
  saveas(gcf,[PATHout,'mapa_topografia'],'png')
%   cd ..
end