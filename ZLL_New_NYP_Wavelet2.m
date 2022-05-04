%%%%%%%%%%%Wavelet for space%%%%WAVETEST Example Matlab script%%%%%
%%%%%%%%%Zi-Liang Li%%%%%2017Year 1 October%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%1-MATLAB读取NC文件%%%%%%%%%%%%%%%%%%%%
% ncread
%1.1 功能
%从NetCDF格式的数据源的变量中读取数据（Read data from variable in NetCDF data source）
%1.2 语法结构
%vardata = ncread(source,varname)
%vardata = ncread(source,varname,start,count,stride)
%1.3 描述
%1.3.1 vardata = ncread(source,varname)
%从数据源中读取变量名为varname的变量
%1.3.2 vardata = ncread((source,varname,start,count,stride)
%（1）start
%varname所指定变量的每一维的开始读取的位置
%（2）count
%从start指定的开始位置算起，一共读取的每一维要素的数目
%（3）stride
%从start开始，每一维读取的数目为count时，每一维的读取的步长
%%%%%%%%%%%%%%1-MATLAB读取NC文件%%%%%%%%%%%%%%%%%%%%
clear;
%%%clear ALL;
clc;
%%%%%%%%%%%%%Figure 1%%%%%%%%%%%%%%%
subplot(2,2,[1,2]);
subplot('position',[0.125 0.605 0.82 0.32]);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%%%%%%%%%%%%%%%%%Figure 2%%%%%%%%%%%%%%%
dataname='./wrfout_d02_2013-07-12_12_00_00.nc';
ncinfo(dataname);
aaa=ncreadatt(dataname,'W','coordinates');
%1%%%%%%%%%%begin ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdisp('./zll_sunwei_F.hdfw2.nc','w');  %获取所读取nc文件的基本信息%%%%%%%%%%%
ncdisp(dataname,'W');

ncid = netcdf.open(dataname, 'NC_NOWRITE');
latd = ncread (dataname,'XLAT');
lond = ncread (dataname,'XLONG');
timed = ncread(dataname,'XTIME');
%lat2=lat(:,2,16);
%squeeze(lon)
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
%distance=ncread(dataname,'XLONG');
distance=lond(1:192,1);%start=65, count=116 thus end=65+116-1=180
%ww=ncread(dataname,'W');
%load("ww.txt");
%%ww = sin(5*distance)+sin(20*distance);
ww=ncread(dataname,'T',[1,95,42,15],[192,1,1,1],[1,1,1,1]);
disp('ww');

%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n=length(distance);
n1=length(ww);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(ww)^2;
ww=(ww - mean(ww))/sqrt(variance);
%w(:)=w(:)- mean(w(:))/sqrt(variance);
%variance = std(w(65:200,126,25))^2;
%w(65:200,126,25)=w(65:200,126,25)- mean(w(65:200,126,25))/sqrt(variance);
%w(:,z,i)=w(:,z,i) - mean(w(:,z,i))/sqrt(variance);
%2%%%%%%%%%%end variance%%%%%%%%%%%%%%%%%%%%%%
%3%%%%begin wavelet transformation%%%%%%%%%%%%%%%%%%%%%%
mother = 'Morlet';  %MOTHER = the mother wavelet function%
                             %The choices are 'MORLET', 'PAUL', or 'DOG'%
dstride         = 1;                            
%dstride         = 1;    %dstride=distance(2)-distance(1);%distance = distance ; %%%
                             %sampling rate DT=dstride%%%%%%%%%
pad       = 0;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = 2*dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wave,wavelength,scale1,coi1] = wavelet(ww(:),dstride,pad,dj,s0,j1,mother);
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
%设置绘图的大小,不需要到word里再调整大小(7cm)%
%set(gca,'Position',[.13 .17 .80 .74]); 
%get hanlde to current axis返回当前图形的当前坐标轴的句柄%
%设置xy轴在图片中占的比例，可能需要自己微调%%
%%%%%Figure resize%%%%%%%%%%%
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
%将字体大小改为8号字，在小图里很清晰%%%
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%将线宽改为2%%'Fontname', 'Times newman'%%
%%%%%%%%%%PLOT%%%%%%%%%%%%
%levels = (0.025:0.25:5.5); 
levels = (-65:0.70:65);  %for color figure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v = [0,0.5,1.0,1.5];
Xtick_pos=1:5:192; %确定Label显示的位置
XtickLabel={'115°E','118°E','121°E','124°E','127°E','130°E','133°E'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yticks=0:100:490;
YtickLabel = 2.^(fix(log2(min(wavelength))):fix(log2(max(wavelength))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ytick_pos=1:10:159; %确定Label显示的位置;Yticks=0:5:159;
%and modify ylim=[0,60];
%%%Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
%[c11,h11] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2); 
[c11,h11] = contourf(distance,log2(wavelength),real_wave,levels); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
%%%title('(a) Latitude:23.76°N Height:30km  Wavelet Power');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('(a)Real Wave');
%%%set(title('(a) Latitude:23.76°N Height:30km WS'),'FontName','Times New Roman','FontSize',18,'Color','b')
set(title('(a)Real Wave'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%%ylim = [0,fix(max(wavelength))];
%%ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%   
set(gca,'YLim',log2([min(wavelength),max(wavelength)]), ...
    'YDir','default', ...
'YTick',log2(YtickLabel(:)), ...
'YTickLabel',log2(YtickLabel));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorDepth = 1000;
colormap(jet(colorDepth));
colorbar;
colorbar('FontSize',18,'Color','b');
%colorbar('location','southoutside');
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
%levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
%[c12,h12] = contour(distance,wavelength,real_wave,levels,'b-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Wavelength%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%绘制小波功率谱等值线图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sst = ww(:); 
dlmwrite('ww.txt',sst,'precision','%5.3f');
%------------------------------------------------------ Computation
% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
%variance =std(sst)^2;
%sst = (sst -mean(sst))/sqrt(variance) ;
n = length(sst);
dt = 1;
time = (0:length(sst)-1)*dt; 
%%%xlim = [0,1920];  %% plotting range
pad = 0;     
dj = 1/12; %%%%%0.125;   
s0 = 2*dt;  
j1 = 6/dj;
%%%j1 = (log2(n*dt/s0))/dj;
lag1 = 0.72;
mother = 'Morlet';
[wave,period,scale,coi]= wavelet(ww(:),dt,pad,dj,s0,j1,mother);
power =(abs(wave)).^2 ;  
[signif,fft_theor]= wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 =(signif')*(ones(1,n));
sig95 = power ./sig95;   
global_ws =variance*(sum(power')/n);
dof = n - scale;
global_signif =wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);
%%绘图      
subplot(2,2,3);
%subplot('position',[0.1 0.15 0.650 0.35]);

set(findobj('FontSize',18),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---绘制小波能量谱
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16];
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
contourf(time,log2(period),log2(power),log2(levels));
%%%[c11,h11] = contourf(distance,wavelength,real_wave,levels); 
%%%%xlabel('年份')
%ylabel('wavelength(km)')
%%%%%title('b) 小波能量谱')
%%%set(gca,'XLim',xlim(:))
%set(gca,'YLim',log2([min(period),max(period)]),...
%         'YDir','default', ...
%         'YTick',log2(Yticks(:)), ...
%         'YTickLabel',Yticks);
imagesc(time,log2(period),log2(power));
%%%set(gca,'XLim',xlim(:))
%%%distance=lond(1:192,1);%start=65, count=116 thus end=65+116-1=180
%%%xlim = [distance(1),distance(length(distance))];
%%%Xtick_pos=1:5:192; %确定Label显示的位置
%%%XtickLabel={'115°E','118°E','121°E','124°E','127°E','130°E','133°E'};
%%%set(gca,'XLim',xlim(:), ...
%%%    'XTick',Xtick_pos(:), ...
%%%    'XTickLabel', XtickLabel);
%%%%set(gca,'YLim',[0,1.25*max(global_ws)],'FontSize',18) 
 set(gca,'YLim',log2([min(period),max(period)]),...
          'YDir','default', ...
         'YTick',log2(Yticks(:)), ...
         'YTickLabel',log2(Yticks));
set(title('b) WS(Morlet)'),'FontName','Times New Roman','FontSize',18,'Color','b')
colorDepth = 1000;
colormap(jet(colorDepth));

xlabel('Distance (km)')
ylabel('wavelength(km)')
title('b)WS(Morlet)')

set(gca,'FontName','Times New Roman','FontSize',18);
%95% singificance contour,levelsat -99(fake)and 1(95% signif)
set(findobj('FontSize',18),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
hold on
contour(time,log2(period),sig95,[-99,1],'k','linewidth',2);
hold on
plot(time,log2(coi),'k','linewidth',2,'Color','b')
%%%%colorbar;
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%绘制全域能谱
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
set(gca,'YLim',[min(log2(period)),max(log2(period))],...
         'YDir','default', ...
         'YTick',log2(Yticks(:)), ...
         'YTickLabel',log2(Yticks));
set(gca,'XLim',[0,1.25*max(global_ws)]);
set(title('c) Global WS'),'FontName','Times New Roman','FontSize',18,'Color','b')
colorDepth = 1000;
colormap(jet(colorDepth));
set(gca,'FontName','Times New Roman','FontSize',18);
set(findobj('FontSize',18),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end of code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
