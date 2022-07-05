% @Mail : chsegura@udec.cl
% @Other : Main source from tskew.m by Steven K. Krueger 

% @Last Modified by:   Christian Segura
% @Last Modified time: 2019-09-05 
% Modificado por MA Troncoso Villar 2022-07-4

function skewt_wrfout(dirdata,filename,xx,yy,t1,inicio)

% posicion grilla a evaluar
% xx=[-71]; % longitud
% yy=[-36]; % latitud

% fecha a evaluar
% t1={'03-Jan-2020 12:00:00'}

%directorio 
% dir='/home/chris/Documentos/Toolbox_WRF_matlab';
% dir='/home/matlab/WRF/Toolbox_WRF_matlab';

% addpath([dir]);
% addpath([dir,'/funciones']);
% addpath([dir,'/scripts']);
% addpath([dir,'/funciones/m_map']);
% addpath(['/home/matlab/croco_tools/UTILITIES/m_map1.4h']);

% inicio='03-Jan-1979 00:00:00';%Hora de inicio simulacion, por ej. '01-Jan-1979 00:00:00'
PATHin = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\exp006\';
PATHout = 'C:\Users\fredo\OneDrive\Documents\Geof\HPTesis\Plots\';

d01=[dirdata,filename];%,'.nc'];


%% No modificar

XLAT=ncread(d01,'XLAT');%latitud (sur negativo)
XLONG=ncread(d01,'XLONG');%longitud (oeste negativo)

XLAT=double(XLAT);
XLAT=squeeze(XLAT(:,:,1));

XLONG=double(XLONG);
XLONG=squeeze(XLONG(:,:,1));

lat=XLAT(1,:)';
lon=XLONG(:,1);

x=lon; % x vector arbitrario
y=lat; % x vector arbitrario

clear minimox positionx minimoy positiony mas_cercanox mas_cercanoy
for i=1:length(xx)
[minimox(i),positionx(i)]=min(abs(x-xx(i)));
[minimoy(i),positiony(i)]=min(abs(y-yy(i)));
mas_cercanox(i)=x(positionx(i));
mas_cercanoy(i)=y(positiony(i));
end
if mas_cercanoy<0
    NorSur = "°S ";
else
    NorSur = "°N ";
end
if mas_cercanox<0
    EstOst = "°O, ";
else
    EstOst = "°E, ";
end
coord_final = " "+abs(fix(100*mas_cercanoy)/100)+NorSur+abs(fix(100*mas_cercanox)/100)+EstOst;

clear lat lon mas_cercanox mas_cercanoy minimox minimoy x XLAT XLONG xx y yy 

%% Recoleccion data 

myVar='XTIME';% minutos desde iniciada la simulación
data=ncread(d01,myVar);
out1 = minutes(diff(datetime([inicio;t1])));

idx=(floor(data)>=out1);

auxnum=find(idx==1);
auxs=auxnum(1);
ts_hour=double(data(auxs)./60);

clear auxnum idx out1 

% myVar='PH';% 
% data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
% PH=double(data);
% PH=squeeze(PH);
% 
% myVar='PHB';% 
% data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
% PHB=double(data);
% PHB=squeeze(PHB);

myVar='P';% 
data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
P=double(data);
P=squeeze(P);

myVar='PB';% 
data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
PB=double(data);
PB=squeeze(PB);

PT=(P+PB);% presion total en Pa

myVar='T';% perturbacion temperatura potencial
data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
T=double(data)+300;% (theta+t0)
T=squeeze(T);% temperatura potencial

TEMP=T.*((PT./100000).^(2/7));% temperatura kelvin
TEMPC=TEMP-273.15;% temperatura celcius

myVar='QVAPOR';% 
data=ncread(d01,myVar,[positionx positiony 1 auxs ],[1 1 Inf 1]);
QVAPOR=double(data);
QVAPOR=squeeze(QVAPOR);

es = 6.112.*exp(17.67.*(TEMP-273.15)./(TEMP-29.65));
e = QVAPOR .* PT./100./(QVAPOR + 0.622);
RH = e./es;% humedad relativa 0-1

pz=PT./100;
tz=TEMPC;
rhz=RH;
% 
uwind=ncread([PATHin,filename],'U',[positionx positiony,1,auxs],[1,1,Inf,1]);
uwind = squeeze(double(uwind)); % (verticales)
u=uwind;
vwind=ncread([PATHin,filename],'V',[positionx positiony,1,auxs],[1,1,Inf,1]);
vwind = squeeze(double(vwind)); % (verticales)
v=vwind;
%
%%
%Now create wind barbs on the right
%Based on this file:
%  MFILE:   windbarbm.m
%  MATLAB:  7.8.0 (R2009a)
%  VERSION: 1.3 (28 November 2011)
%  AUTHOR:  Nick Siler
%  CONTACT: siler@atmos.washington.edu
% if ~isempty(u)
    %This is important only if wind barbs are displayed
%     axis square
    ctp=100; % presión más baja (tope del plot)
    u_barb=u(pz>=ctp);
    v_barb=v(pz>=ctp);
%     u = u.*1.944; %pasarlo a nudos????
%     v = v.*1.944;
    pz_barb=pz(pz>=ctp);
    
    offset=-48;  %-45; %Minimum chart T
    nj=50-offset+5; %Range of T

    xloc = offset+nj+3; %change this to change the position of the wind barb axis
    
    cbp=1010; %presión más alta (superficie)
    xx=ones(size(u_barb))*xloc;
    
T0=273.15;
R=287;
Lv=2.5e6;
H=10;
%
    yy=H*log(cbp./pz_barb);
    
        umag = (u_barb.^2+v_barb.^2).^0.5; %wind speed
        %find theta; add pi to atan(v/u) when u<0
        dummy = (u_barb<0)*pi;
        theta_dir = atan(v_barb./u_barb)+dummy;
    clear dummy
    [a,b] = size(umag);
    
    %create 18 logical matrices for 18 possible barbs. Non-zero when the barb
    %is called for at that gridpoint.
    g=cell(1,18);
    g{1} = umag > 7.5 & umag <= 47.5;
    g{2} = umag > 17.5 & umag <= 47.5;
    g{3} = umag > 27.5;
    g{4} = (umag > 37.5 & umag <= 47.5) | (umag > 57.5 & umag <= 97.5);
    g{5} = umag > 67.5;
    g{6} = (umag > 77.5 & umag < 97.5) | umag > 107.5;
    g{7} = umag > 87.5 & umag < 97.5 | umag > 117.5;
    g{8} = umag > 127.5;
    g{9} = (umag > 2.5 & umag <= 7.5) | (umag > 12.5 & umag <= 17.5);
    g{10} = umag > 22.5 & umag <= 27.5;
    g{11} = (umag > 32.5 & umag <= 37.5) | (umag > 52.5 & umag <= 57.5);
    g{12} = (umag > 42.5 & umag <= 47.5) | (umag > 62.5 & umag <= 67.5);
    g{13} = (umag > 72.5 & umag <= 77.5) | (umag > 102.5 & umag <= 107.5);
    g{14} = (umag > 82.5 & umag <= 87.5) | (umag > 112.5 & umag <= 117.5);
    g{15} = (umag > 92.5 & umag <= 97.5) | (umag > 122.5 & umag <= 127.5);
    g{16} = umag > 47.5;
    g{17} = umag > 97.5;
    g{18} = true(a,b);
    
    %position of each barb relative to grid point: [x0 y0; x1 y1]
    c{1} = [-1 0;-1.125 .325];
    c{2} = [-.875 0; -1 .325];
    c{3} = [-.75 0; -.875 .325];
    c{4} = [-.625 0; -.75 .325];
    c{5} = [-.5 0; -.625 .325];
    c{6} = [-.375 0; -.5 .325];
    c{7} = [-.25 0; -.375 .325];
    c{8} = [-.125 0; -.25 .325];
    c{9} = [-.875 0; -.9375 .1625];
    c{10} = [-.75 0; -.8125 .1625];
    c{11} = [-.625 0; -.6875 .1625];
    c{12} = [-.5 0; -.5625 .1625];
    c{13} = [-.3750 0; -.4375 .1625];
    c{14} = [-.25 0; -.3125 .1625];
    c{15} = [-.125 0; -.1875 .1625];
    c{16} = [-1 0; -.875 .325];
    c{17} = [-.75 0; -.625 .325];
    c{18} = [0 0; -1 0];
    
    for nn=1:18
        c{nn}=c{nn}*2;
    end
    
    scale2x=3;
    scale2y=scale2x*diff(H*log(cbp./[cbp ctp]))/(nj-5);

%% definiciones skewT
ez=6.112.*exp(17.67.*tz./(243.5+tz));
qz=rhz.*0.622.*ez./(pz-ez);
chi=log(pz.*qz./(6.112.*(0.622+qz)));
tdz=243.5.*chi./(17.67-chi);
%
p=[1050:-25:100];
pplot=transpose(p);
t0=[-48:2:50];
[ps1,ps2]=size(p);
ps=max(ps1,ps2);
[ts1,ts2]=size(t0);
ts=max(ts1,ts2);
for i=1:ts,
   for j=1:ps,
      tem(i,j)=t0(i)+30.*log(0.001.*p(j));
      thet(i,j)=(273.15+tem(i,j)).*(1000./p(j)).^.287;
      es=6.112.*exp(17.67.*tem(i,j)./(243.5+tem(i,j)));
      q(i,j)=622.*es./(p(j)-es);
      thetaea(i,j)=thet(i,j).*exp(2.5.*q(i,j)./(tem(i,j)+273.15));
   end
end
%%
tname="SkewT" + coord_final + datestr(datetime(t1),'yyyy/mm/dd HH:MM');% titulo de la figura
figure
p=transpose(p);
t0=transpose(t0);
temp=transpose(tem);
theta=transpose(thet);
thetae=transpose(thetaea);
qs=transpose(sqrt(q));
h=contour(t0,pplot,temp,16,'k'); %negras isotermas
hold on
h=contour(t0,pplot,theta,24,'b'); %azules adiabatica seca
% levels = contourdata(h);
% levels = extractfield(levels,'level');
% clabel(h,fix(100.*levels(1:3:end))/100, 'labelspacing', 1500,'Color','b')
h=contour(t0,pplot,thetae,24,'r'); %rojas adiabatica húmeda
% levels = contourdata(h);
% levels = extractfield(levels,'level');
% clabel(h,fix(10000.*levels(2:3:end))/10000, 'labelspacing', 1500,'Color','r')
[h,cc]=contour(t0,pplot,qs,[0.1 0.5 1 1.5 2 3 5 7 9 12],'color','#77AC30','showtext','on'); %el de arriba pero con labels fijas
clabel(h,cc,'color',[0.4660 0.6740 0.1880])% levels = contourdata(h);
% h=contour(t0,pplot,qs,24,'color','#77AC30'); %verdes razon de mezcla del vapor de agua
% levels = extractfield(levels,'level');
% clabel(h,levels(1:2:end), 'labelspacing', 1500,'Color','#77AC30')
%tsound=30.+43.5.*log(0.001.*p);
%tsoundm=tsound-30.*log(0.001.*p);
tzm=tz-30.*log(0.001.*pz);
tdzm=tdz-30.*log(0.001.*pz);
h=plot(tzm,pz,'r',tdzm,pz,'b--'); %temps

    for nn = 1:18
        dummy = reshape(g{nn},1,a*b);
        count = sum(dummy); % number of barbs to draw
        if count == 0
            continue
        end
        
        %rotation operations
        x1 = c{nn}(1,1)*cos(theta_dir)-c{nn}(1,2)*sin(theta_dir);
        y1 = c{nn}(1,1)*sin(theta_dir)+c{nn}(1,2)*cos(theta_dir);
        x2 = c{nn}(2,1)*cos(theta_dir)-c{nn}(2,2)*sin(theta_dir);
        y2 = c{nn}(2,1)*sin(theta_dir)+c{nn}(2,2)*cos(theta_dir);
        
        x1 = x1*scale2x+xx;
        x2 = x2*scale2x+xx;
        
        y1 = y1*scale2y+yy;
        y2 = y2*scale2y+yy;
        
        x = [reshape(x1(dummy),1,count);reshape(x2(dummy),1,count)];
        y = [reshape(y1(dummy),1,count);reshape(y2(dummy),1,count)];
        y = cbp*exp(-y/H);
        
       line(x,y,'color','k','linestyle','-','linewidth',1,'clipping','off','marker','none')
        if nn==18
            plot(x(1,:),y(1,:),'o','color','k','MarkerSize',2,'clipping','off')
        end
    end
    hold on
    plot([xloc xloc],[ctp cbp],'color','k','linestyle','-','linewidth',1,'clipping','off')
hold on
set(gca,'ytick',[1000:100:100])
set(gca,'yscale','log','ydir','reverse')
set(gca,'fontweight','bold')
set(gca,'ytick',[100:100:1000])
set(gca,'ygrid','on')
set(h,'linewidth',2)
% hold off
xlabel('Temperatura [°C]','fontweight','bold')
ylabel('Presión [hPa]','fontweight','bold')
legend(h,'T parcela','T p. de rocío','location','southwest')
ylim([100 1000])
title(tname)
% % % 
  plotname=['skewt_',datestr(datetime(t1),'yymmdd_HHMM')];% nombre output
% % saveas(h,[pwd,'/../figuras/',plotname],'png');
  saveas(gcf,[PATHout,plotname],'png')
% close all

return
