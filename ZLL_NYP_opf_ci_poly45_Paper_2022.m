%%%%%%%%%%%Wavelet for space%%%%WAVETEST Example Matlab script%%%%%
%%%%%%%%%Zi-Liang Li%%%%%2017Year 1 October%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%1-MATLAB��ȡNC�ļ�%%%%%%%%%%%%%%%%%%%%
% ncread
%1.1 ����
%��NetCDF��ʽ������Դ�ı����ж�ȡ���ݣ�Read data from variable in NetCDF data source��
%1.2 �﷨�ṹ
%vardata = ncread(source,varname)
%vardata = ncread(source,varname,start,count,stride)
%1.3 ����
%1.3.1 vardata = ncread(source,varname)
%������Դ�ж�ȡ������Ϊvarname�ı���
%1.3.2 vardata = ncread((source,varname,start,count,stride)
%��1��start
%varname��ָ��������ÿһά�Ŀ�ʼ��ȡ��λ��
%��2��count
%��startָ���Ŀ�ʼλ������һ����ȡ��ÿһάҪ�ص���Ŀ
%��3��stride
%��start��ʼ��ÿһά��ȡ����ĿΪcountʱ��ÿһά�Ķ�ȡ�Ĳ���
%%%%%%%%%%%%%%1-MATLAB��ȡNC�ļ�%%%%%%%%%%%%%%%%%%%%
clear;
%%%clear ALL;
clc;
%%%%%%%%%%%%%Figure 1%%%%%%%%%%%%%%%
%subplot(2,2,1);
%subplot('position',[0.125 0.62 0.78 0.32]);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%%%%%%%%%%%%%%%%%Figure 2%%%%%%%%%%%%%%%
%%load 'dataopf.dat';
%%load 'datapoly4.dat';
load 'datacress.dat';
%%wwa=dataopf;
%%wwa=datapoly4;
wwa=datacress;
%%ncinfo(dataname);
%%aaa=ncreadatt(dataname,'W','coordinates');
%1%%%%%%%%%%begin ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdisp('./zll_sunwei_F.hdfw2.nc','w');  %��ȡ����ȡnc�ļ��Ļ�����Ϣ%%%%%%%%%%%
%%ncdisp(dataname,'W');
%%ncid = netcdf.open(dataname, 'NC_NOWRITE');
%%latd = ncread (wwa,'XLAT');
%%lond = ncread (wwa,'XLONG');
%%timed = ncread(dataname,'XTIME');
%lat2=lat(:,2,16);
%squeeze(lon)
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
%distance=ncread(dataname,'XLONG');
%%distance2=lond(1:192,80,10);
%%distance=lond(1:192,1);%start=65, count=116 thus end=65+116-1=180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ww=ncread(dataname,'W');
%%wwa=load("mountian_w.dat");
%%%wwa = load('ww.txt');
%wwa = sin(5*distance2)+sin(20*distance2);
%%wwa=ncread(dataname,'W',[1,70,42,15],[192,1,1,1],[1,1,1,1]);
%%wwa=ncread(dataname,'W',[1,70,42,15],[192,1,1,1],[1,1,1,1]);
%%%[var1,~] = Get_InterpVar(dataname,'W',26,16,0,'z');
%%wwname=ncread(dataname,'W',[1,1,42,15],[192,159,1,1],[1,1,1,1]);
%wwname=wname-mean(wname);
%%disp('wwa');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%W-velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%set(gcf,'Position',[100 100 260 220]);
%subplot(3,2,1);
%%TightPlot.ColumeNumber = 3;     % ��ͼ����
%%TightPlot.RowNumber = 2;    % ��ͼ����
%%TightPlot.GapW = 0.08;  % ��ͼ֮������Ҽ��
%%TightPlot.GapH = 0.1;   % ��ͼ֮������¼��
%%TightPlot.MarginsLower = 0.1;   % ��ͼ��ͼƬ�Ϸ��ļ��
%%TightPlot.MarginsUpper = 0.03;  % ��ͼ��ͼƬ�·��ļ��
%%TightPlot.MarginsLeft = 0.06;   % ��ͼ��ͼƬ�󷽵ļ��
%%TightPlot.MarginsRight = 0.01;  % ��ͼ��ͼƬ�ҷ��ļ��
%%pp = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
%%    [TightPlot.GapH TightPlot.GapW],...
%%    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
%%    [TightPlot.MarginsLeft TightPlot.MarginsRight]);
%%axes(pp);
%%set(pp,'XTickLabel','');    %%Ĩȥ��ͼ1-4�ĺ�����ֵ
%%set(pp,'YTickLabel','') %%Ĩȥ��ͼ2��4��6��������ֵ
%%%set(gca,'position', [0.05 0.15 0.3 0.58]);
%%lond2=ncread(dataname,'XLONG');latd2=ncread(dataname,'XLAT');
%%lon2=lond2(:,:,16);lat2=latd2(:,:,16);
%%latmin = 17;
%%latmax = 30;
%latint = 10;
%lonint = 10;
%%lonmin = 115;
%%lonmax =132;
%%yt = [17 23 30];
%%xt = [115 123 132];
%%xtl=['115��E';'123��E';'132��E'];
%%ytl=['17��N';'23��N';'30��N'];
%%left = [.06,.545];
%%bottom = [.59,.17];
%%width = 0.615;
%%height = .325;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
%%m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',12,...
%%'tickdir','out','Color','b','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
%%hold on;
%%m_coast('color',[0.7 0.7 0.7]);hold on;
%%dlevels = [-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4] ;
%%dlevels = [-5,-2,-0.5,-0.1,0,0.1,0.5,2,5];
%%wname=ncread(dataname,'W',[1,1,42,15],[192,159,1,1],[1,1,1,1]);
%%wnamenew=wname-mean(wname);
%%wnamez=makecolor(wnamenew,dlevels);
%%m_contourf(lon2,lat2,wnamez,'linestyle','none');
%%hold on;%m_contourf(long1,lati,real(tbb));
%%m_coast;
%%dlevel = [-5,0,5];
%%h1= colorbar('FontSize',12,'Color','b');% 'FontWeight', 'bold',
%%set(h1,'Ticks',[0,4,8],'TickLabels',dlevel) ;
%set(get(h1,'ylabel'),'string','w/m.s^-^1','fontname','Times New Roman',...
%    'fontsize',16);%,[a b]Ϊ��ͼ�������µ�����ꡣc��d�ֱ�Ϊ��ͼ����Ŀ�͸ߡ� 
%colorbar;% 'FontWeight', 'bold',
%colorbar('FontSize',18,'Color','b');
%hc = colorbar;colormap([1 0 0;0 1 0])
%set(hc,'YTick',dlevels);
%set(hc,'YTickLabel',[min(min(wname)) 0 max(max(wname))]);
%colorbar('Ticks',[-5,-2,1,4,7],...
%         'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
%%set(gcf,'Position',[100 100 260 220]);
%%title('(a)Vertical Velocity');
%%%set(title('(a) Latitude:23.76��N Height:30km WS'),'FontName','Times New Roman','FontSize',18,'Color','b')
%%set(title('a)Vertical Velocity'),'FontName','Times New Roman','FontSize',12,'Color','b')
%%grid on
%%hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%subplot(3,2,2);
%%lond2=ncread(dataname,'XLONG');latd2=ncread(dataname,'XLAT');
%%lon2=lond2(:,:,16);lat2=latd2(:,:,16);
%%latmin = 17;
%%latmax = 30;
%%latint = 10;
%%lonint = 10;
%%lonmin = 115;
%%lonmax =132;
%%yt = [17 23 30];
%%xt = [115 123 132];
%%xtl=['115��E';'123��E';'132��E'];
%%ytl=['17��N';'23��N';'30��N'];
%left = [.06,.545];
%bottom = [.59,.17];
%width = .415;
%height = .325;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
%%m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',12,...
%%'tickdir','out','Color','b','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
%%hold on;
%%m_coast('color',[0.7 0.7 0.7]);hold on;
%%dlevels = [-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4] ;
%%w2name=ncread(dataname,'T',[1,1,42,15],[192,159,1,1],[1,1,1,1]);
%%w2namet=w2name-mean(w2name);
%%dlevels = [-10,-8,-6,-4,-2,0,2,4,6,8,10];
%%w2namez=makecolor(w2namet,dlevels);
%%m_contourf(lon2,lat2,w2namez,'linestyle','none');
%%hold on;%m_contourf(long1,lati,real(tbb));
%%m_coast;
%%title('(b)Vertical Velocity');
%%%set(title('(a) Latitude:23.76��N Height:30km WS'),'FontName','Times New Roman','FontSize',18,'Color','b')
%%set(title('b)PT Perturbation'),'FontName','Times New Roman','FontSize',12,'Color','b')
%colorbar;% 'FontWeight', 'bold',
%colorbar('FontSize',18,'Color','b');
%%dlevel = [-5,0,5];
%%cbwt = .49;
%%h1= colorbar('FontSize',12,'Color','b');% 'FontWeight', 'bold',
%%set(h1,'Ticks',[0,4,8],'TickLabels',dlevel);
%%hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end w-velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count=[136,1,1]��ÿһά��start��ʼ��ȡ������Ŀ%%stride=[1,1,1]���ö�ȡ�Ĳ���%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n=length(distance);
subplot(2,2,[1,2]);
nlength=length(wwa);
dt=1.7184;
timedt = (0:(nlength-1))*dt; 
distance = timedt + 0.0; 
%%distance=lond(1:n1);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variance = std(wwa)^2;
%%%%%%%%%%%%%%���㷽��%%%%%%%%%%%%%%%%%%%%
wwa=(wwa - mean(wwa))/sqrt(variance);
%ww = smoothdata(wwq,'gaussian',10);%����aΪ���ݣ�'gaussian'Ϊ��˹�˲�����
%w(:)=w(:)- mean(w(:))/sqrt(variance);
%variance = std(w(65:200,126,25))^2;
%w(65:200,126,25)=w(65:200,126,25)- mean(w(65:200,126,25))/sqrt(variance);
%w(:,z,i)=w(:,z,i) - mean(w(:,z,i))/sqrt(variance);
%2%%%%%%%%%%end variance%%%%%%%%%%%%%%%%%%%%%%
%3%%%%begin wavelet transformation%%%%%%%%%%%%%%%%%%%%%%
mother = 'Morlet';  %MOTHER = the mother wavelet function%
                             %The choices are 'MORLET', 'PAUL', or 'DOG'%
dstride         = 1.7184;                            
%dstride         = 1;    %dstride=distance(2)-distance(1);%distance = distance ; %%%
                             %sampling rate DT=dstride%%%%%%%%%
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = 2*dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
%%%j1          = 7.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
j1 = fix((log2(nlength*dstride/s0))/dj);
%%%lag1      = 0.72;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wave,wavelength,scale1,coi1] = wavelet(wwa(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WAVE is the WAVELET transform of Y. This is a complex array
%of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%The WAVELET power spectrum is ABS(WAVE)^2.
%Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%
%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set(gcf,'Position', [0.505 0.505 0.495 0.495]);
%���û�ͼ�Ĵ�С,����Ҫ��word���ٵ�����С(7cm)%
%set(gca,'Position',[.13 .17 .80 .74]); 
%get hanlde to current axis���ص�ǰͼ�εĵ�ǰ������ľ��%
%����xy����ͼƬ��ռ�ı�����������Ҫ�Լ�΢��%%
%%%%%Figure resize%%%%%%%%%%%
figure_FontSize=12;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',12),'FontSize',figure_FontSize);
%�������С��Ϊ8���֣���Сͼ�������%%%
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%���߿��Ϊ2%%'Fontname', 'Times newman'%%
%%%%%%%%%%PLOT%%%%%%%%%%%%
%levels = (0.025:0.25:5.5); 
levels = (-50:3:50);  %for color figure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v = [0,0.5,1.0,1.5];
Xtick_pos=1:80:nlength; %ȷ��Label��ʾ��λ��
XtickLabel={'115��E','118��E','121��E','124��E','127��E','130��E','133��E'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yticks=0:100:490;
YtickLabel = 2.^(fix(log2(min(wavelength))):fix(log2(max(wavelength))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ytick_pos=1:10:159; %ȷ��Label��ʾ��λ��;Yticks=0:5:159;
%and modify ylim=[0,60];
%%%Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
%[c11,h11] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2); 
[c11,h11] = contourf(distance,log2(wavelength),real_wave,levels); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
%%%title('(a) Latitude:23.76��N Height:30km  Wavelet Power');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('(a)Real Wave');
%%%set(title('(a) Latitude:23.76��N Height:30km WS'),'FontName','Times New Roman','FontSize',18,'Color','b')
set(title('a)Real Wave'),'FontName','Times New Roman','FontSize',12,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',12,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%%ylim = [0,fix(max(wavelength))];
%%ylim=[2,92];%����������������
%%set(gca,'XLim',xlim(:), ...
%%    'XTick',Xtick_pos(:), ...
%%    'XTickLabel', XtickLabel);
set(gca,'YLim',log2([min(wavelength),max(wavelength)]),...
          'YDir','default', ...
         'YTick',log2(YtickLabel(:)), ...
         'YTickLabel',YtickLabel);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorDepth = 1000;
colormap(jet(colorDepth));
colorbar;
colorbar('FontSize',12,'Color','b');
%colorbar('location','southoutside');
set(gca,'FontName','Times New Roman','FontSize',12);
%���������壥
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
%levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
%[c12,h12] = contour(distance,wavelength,real_wave,levels,'b-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X�����ɫ%
set(gca, 'YColor','blue');        % Y�����ɫ%
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Wavelength%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%����С�������׵�ֵ��ͼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sst = wwa(:); 
%%%dlmwrite('ww.txt',sst,'precision','%5.3f');
%------------------------------------------------------ Computation
% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
%variance =std(sst)^2;
%sst = (sst -mean(sst))/sqrt(variance) ;
n = length(sst);
dt = 1.7184;
time = (0:length(sst)-1)*dt; 
%%%xlim = [0,1920];  %% plotting range
pad = 1;     
dj = 1/12; %%%%%0.125;   
s0 = 2*dt;  
%%j1 = 7.5/dj;
j1 = fix((log2(nlength*dstride/s0))/dj);
%%%j1 = (log2(n*dt/s0))/dj;
lag1 = 0.72;
mother = 'Morlet';
[wave,period,scale,coi]= wavelet(wwa(:),dt,pad,dj,s0,j1,mother);
power =(abs(wave)).^2 ;  
[signif,fft_theor]= wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 =(signif')*(ones(1,n));
sig95 = power ./sig95;   
global_ws =variance*(sum(power')/n);
dof = n - scale;
global_signif =wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);
%%��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot('position',[0.1 0.15 0.650 0.35]);

set(findobj('FontSize',12),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---����С��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16];
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contourf(time,log2(period),log2(power),log2(levels));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%[c11,h11] = contourf(distance,wavelength,real_wave,levels); 
%%%%xlabel('���')
%ylabel('wavelength(km)')
%%%%%title('b) С��������')
%%%set(gca,'XLim',xlim(:))
%set(gca,'YLim',log2([min(period),max(period)]),...
%         'YDir','default', ...
%         'YTick',log2(Yticks(:)), ...
%         'YTickLabel',Yticks);
imagesc(time,log2(period),log2(power));
%%%set(gca,'XLim',xlim(:))
%%%distance=lond(1:192,1);%start=65, count=116 thus end=65+116-1=180
%%%xlim = [distance(1),distance(length(distance))];
%Xtick_pos=1:5:192; %ȷ��Label��ʾ��λ��
%XtickLabel={'115��E','118��E','121��E','124��E','127��E','130��E','133��E'};
%set(gca,'XLim',xlim(:), ...
%    'XTick',Xtick_pos(:), ...
%    'XTickLabel', XtickLabel);
%%%%set(gca,'YLim',[0,1.25*max(global_ws)],'FontSize',18)
%%%%set(gca,'YLim',[0,1.25*max(global_ws)],'FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 set(gca,'YLim',log2([min(period),max(period)]),...
          'YDir','default', ...
         'YTick',log2(Yticks(:)), ...
         'YTickLabel',Yticks);   
  %set(gca,'YLim',log2([min(period),max(period)]), ...
	%'YDir','reverse', ...
	%'YTick',log2(Yticks(:)), ...
	%'YTickLabel',log2(Yticks)));   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(title('b) WPS(Morlet)'),'FontName','Times New Roman','FontSize',12,'Color','b')
colorDepth = 1000;
colormap(jet(colorDepth));
xlabel('Distance (km)')
ylabel('wavelength(km)')
title('b)WPS(Morlet)')
set(gca,'FontName','Times New Roman','FontSize',12);
%95% singificance contour,levelsat -99(fake)and 1(95% signif)
set(findobj('FontSize',12),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
hold on
contour(time,log2(period),sig95,[-99,1],'k','linewidth',2);
hold on
plot(time,log2(coi),'k','linewidth',2,'Color','b')
%%%%colorbar;
set(gca, 'XColor','blue');       % X�����ɫ%
set(gca, 'YColor','blue');        % Y�����ɫ%
hold off
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ȫ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4);
%subplot('position',[0.75 0.15 0.18 0.35]);
%%plot(global_ws,log2(period))
plot(global_ws,log2(period));
hold on
%%plot(global_signif,log2(period),'--')
plot(global_signif,log2(period),'--')
hold off
xlabel('Power(K^2)')
title('Global WS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'YLim',[min(log2(period)),max(log2(period))],...
         'YDir','default', ...
         'YTick',log2(Yticks(:)), ...
         'YTickLabel',Yticks);
%set(gca,'YLim',log2([min(period),max(period)]), ...
%	'YDir','reverse', ...
%	'YTick',log2(Yticks(:)), ...
%	'YTickLabel',Yticks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'XLim',[0,1.25*max(global_ws)]);
set(title('c) Global WPS'),'FontName','Times New Roman','FontSize',12,'Color','b')
colorDepth = 1000;
colormap(jet(colorDepth));
set(gca,'FontName','Times New Roman','FontSize',12);
set(findobj('FontSize',12),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca, 'XColor','blue');       % X�����ɫ%
set(gca, 'YColor','blue');        % Y�����ɫ%
minpower = (min(log2(power)) + max(log2(power)))*0.5;
%colorDepth = 1000;
%colormap(jet(colorDepth));
colorbar;
colorbar('FontSize',12,'Color','b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(gcf,'zll-datacress-wavelet.eps','-depsc2', '-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end of code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
