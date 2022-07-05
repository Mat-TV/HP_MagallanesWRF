%Función para graficar la disposición de dominios de WRFOUT 
  % dentro de la función están las variables PATH{in,out}
%Argumentos de entrada:
    %filename: nombre del archivo
%Argumentos de salida:
    %ninguno: nada
%  2022-07-01 - Mat Troncoso Villar (matro1432@gmail.com)
function wrfout_dominios(filename1,filename2,filename3)
% aseo
% filename1='wrfout_d01_2020-01-01';
% filename2='wrfout_d02_2020-01-02';
% filename3='wrfout_d03_2020-01-03';
PATHin = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\exp006\';
PATHout = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\Plots\';
lat1=ncread([PATHin,filename1],'XLAT');
lat1 = lat1(:,:,1);
lon1=ncread([PATHin,filename1],'XLONG');
lon1 = lon1(:,:,1);
hgt1=ncread([PATHin,filename1],'HGT',[1,1,1],[Inf,Inf,1]);
% 
lat2=ncread([PATHin,filename2],'XLAT');
lat2 = lat2(:,:,1);
lon2=ncread([PATHin,filename2],'XLONG');
lon2 = lon2(:,:,1);
hgt2=ncread([PATHin,filename2],'HGT',[1,1,1],[Inf,Inf,1]);
% 
lat3=ncread([PATHin,filename3],'XLAT');
lat3 = lat3(:,:,1);
lon3=ncread([PATHin,filename3],'XLONG');
lon3 = lon3(:,:,1);
hgt3=ncread([PATHin,filename3],'HGT',[1,1,1],[Inf,Inf,1]);
%
%% 
figure
hold on
m_proj('lambert','lon',[min(min(double(lon1)))+1 max(max(double(lon1)))-1],...
        'lat',[min(min(double(lat1)))+1 max(max(double(lat1)))-1]);
m_pcolor(lon1(:,:),lat1(:,:),hgt1(:,:));
colormap([m_colmap('blues',5);m_colmap('gland',56)]);
m_pcolor(lon2(:,:),lat2(:,:),hgt2(:,:));
colormap([m_colmap('blues',1);m_colmap('gland',56)]);
m_pcolor(lon3(:,:),lat3(:,:),hgt3(:,:));
colormap([m_colmap('blues',2);m_colmap('gland',56)]);
ax1=m_contfbar(.97,[.5 .9],hgt1(:,:), ...
    [64],'edgecolor','none','endpiece','no');
m_plot([min(min(double(lon2))) min(min(double(lon2)))],[min(min(double(lat2))) max(max(double(lat2)))],'--black', 'linewidth',1.8)
m_plot([max(max(double(lon2))) max(max(double(lon2)))],[min(min(double(lat2))) max(max(double(lat2)))],'--black', 'linewidth',1.8)
m_plot([min(min(double(lon2))) max(max(double(lon2)))],[max(max(double(lat2))) max(max(double(lat2)))],'--black', 'linewidth',1.8)
m_plot([min(min(double(lon2))) max(max(double(lon2)))],[min(min(double(lat2))) min(min(double(lat2)))],'--black', 'linewidth',1.8)
m_text(max(max(double(lon2)))-0.5,max(max(double(lat2)))+0.5,'d02','fontsize',8)
%
m_plot([min(min(double(lon3))) min(min(double(lon3)))],[min(min(double(lat3))) max(max(double(lat3)))],'--black', 'linewidth',1.2)
m_plot([max(max(double(lon3))) max(max(double(lon3)))],[min(min(double(lat3))) max(max(double(lat3)))],'--black', 'linewidth',1.2)
m_plot([min(min(double(lon3))) max(max(double(lon3)))],[max(max(double(lat3))) max(max(double(lat3)))],'--black', 'linewidth',1.2)
m_plot([min(min(double(lon3))) max(max(double(lon3)))],[min(min(double(lat3))) min(min(double(lat3)))],'--black', 'linewidth',1.2)
m_text(max(max(double(lon3)))-0.5,max(max(double(lat3)))+0.5,'d03','fontsize',8)
%
title(ax1,{'Altitud [m]',''}); % Move up by inserting a blank line
m_grid('box','fancy','tickdir','in');
title('Topografía y dominios.')
xlabel('Longitud','fontsize',12)
ylabel('Latitud','fontsize',12)
hold off
  saveas(gcf,[PATHout,'dominios'],'png')
end