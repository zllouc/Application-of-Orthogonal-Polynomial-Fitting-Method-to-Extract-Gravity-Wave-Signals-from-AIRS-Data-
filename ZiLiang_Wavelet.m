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
%subplot(3,2,1);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%%%%%%%%%%%%%Figure 1%%%%%%%%%%%%%%%

%1%%%%%%%%%%begin ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdisp('./zll_sunwei_F.hdfw2.nc','w');  %获取所读取nc文件的基本信息%%%%%%%%%%%
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
%distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 65, 136, 1);
%w=ncread('./zll_sunwei_F.hdfw2.nc','w',[65,129,15],[136,1,1],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
distance=lond(1:192);
%n=length(distance);
w = sin(5*distance)+sin(20*distance);
n1=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);

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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;        
[wave,wavelength,scale1,coi1] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
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
levels = (-6.5:0.70:6.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=65:33:201; %确定Label显示的位置
XtickLabel={ '113°E','114°E', ...
    '115°E', '116°E', '117°E'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:10:92; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=6:30:276; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
%[c11,h11] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2); 
[c11,h11] = contourf(distance,wavelength,real_wave,levels); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('(a) Latitude:26°N Height:7.5km  Wavelet Power');
set(title('(a) Latitude:26°N Height:7.5km Wavelet Power'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
colorDepth = 1000;
colormap(jet(colorDepth));
colorbar;
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

%%%%%%%%%%%%%%%%%Figure 2%%%%%%%%%%%%%%%
subplot(3,2,2);
%set(gca, 'Units', 'normalized', 'Position', [0.505 0.505 0.495 0.495]);
%%%%%%%%%%%%%%%%%Figure 2%%%%%%%%%%%%%%%

%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 65, 136, 1);
w=ncread('./zll_sunwei_F.hdfw2.nc','w',[65,129,20],[136,1,1],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n=length(distance);
n2=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);
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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;        
[wave,wavelength,scale2,coi2] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
                            %WAVE is the WAVELET transform of Y. This is a complex array
                            %of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
                            %ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
                            %The WAVELET power spectrum is ABS(WAVE)^2.
                            %Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%

%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'Position', [400,100,300,300]);
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
levels = (-6.5:0.70:6.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=65:33:201; %确定Label显示的位置
XtickLabel={ '113°E','114°E', ...
    '115°E', '116°E', '117°E'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:10:92; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=6:30:276; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
[c21,h21] = contourf(distance,wavelength,real_wave,levels); 
%[c21,h21] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2);
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('(b) Latitude:26°N Height:10 km Wavelet Power');
set(title('(b) Latitude:26°N Height:10 km Wavelet Power'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
colorDepth = 1000;
colormap(jet(colorDepth));
colorbar;
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
%levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
%[c22,h22] = contour(distance,wavelength,real_wave,levels,'b-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off
%colormap(gca,'jet')
%colorbar
%caxis([-1,1]);

%%%%%%%%%%%%Figure 3%%%%%%%%%%%%%
subplot(3,2,3);
%%%%%%%%%%%%Figure 3%%%%%%%%%%%%%
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 65, 136, 1);
w=ncread('./zll_sunwei_F.hdfw2.nc','w',[65,129,25],[136,1,1],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n=length(distance);
n3=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);
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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;        
[wave,wavelength,scale3,coi3] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
                            %WAVE is the WAVELET transform of Y. This is a complex array
                            %of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
                            %ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
                            %The WAVELET power spectrum is ABS(WAVE)^2.
                            %Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%

%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'Position', [400,100,300,300]);
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
levels = (-6.5:0.70:6.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=65:33:201; %确定Label显示的位置
XtickLabel={ '113°E','114°E', ...
    '115°E', '116°E', '117°E'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:10:92; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=6:30:276; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
%[c31,h31] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2); 
[c31,h31] = contourf(distance,wavelength,real_wave,levels); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('(c) Latitude:26°N Height:12.5km Wavelet Power');
set(title('(c) Latitude:26°N Height:12.5km Wavelet Power'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
colorDepth = 1000;
colormap(jet(colorDepth));
h=colorbar;
%set(get(h,'Title'),'string','m');
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
%levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
%[c32,h32] = contour(distance,wavelength,real_wave,levels,'b-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off

%%%%%%%%%%%Figure 4%%%%%%%%%%%%%%%%%%%%
subplot(3,2,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 1, 36, 1);
w=ncread('./zll_sunwei_F.hdfw2.nc','w',[159,129,1],[1,1,36],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n=length(distance);
n4=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);
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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;        
[wave,wavelength,scale4,coi4] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
                            %WAVE is the WAVELET transform of Y. This is a complex array
                            %of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
                            %ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
                            %The WAVELET power spectrum is ABS(WAVE)^2.
                            %Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%

%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'Position', [400,100,300,300]);
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
levels = (-5.5:0.65:5.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=1:6:36; %确定Label显示的位置
XtickLabel={ '3','6', ...
    '9', '12', '15','Height(km)'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:6:36; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=1:3:18; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
[c41,h41] = contourf(distance,wavelength,real_wave,levels); 
%[c41,h41] = contourf(distance,wavelength,real_wave,levels,'r--','linewidth',1.2); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('高度:12km');
set(title('高度:12km'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel('')
%set(xlabel('Height(km)'),'FontSize',18,'Color','b')
ylabel('Wavelength(km)')
set(ylabel('Wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,36];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
colorDepth = 1000;
colormap(jet(colorDepth));
colorbar;
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
%levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
%[c42,h42] = contour(distance,wavelength,real_wave,levels,'b-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off

%%%%%%%%%%%%Figure 5%%%%%%%%%%%%%
subplot(3,2,5);
%%%%%%%%%%%%Figure 5%%%%%%%%%%%%%
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
%distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
%w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 65, 136, 1);
w=ncread('./zll_sunwei_F.hdfw2.nc','w',[65,126,7],[136,1,1],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n=length(distance);
n5=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);
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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
%lag1      = 0.72;        
[wave,wavelength,scale5,coi5] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
                            %WAVE is the WAVELET transform of Y. This is a complex array
                            %of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
                            %ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
                            %The WAVELET power spectrum is ABS(WAVE)^2.
                            %Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%

%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'Position', [400,100,300,300]);
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
levels = (0.025:0.25:5.5); 
%levels = (-5.5:0.065:5.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=65:33:201; %确定Label显示的位置
XtickLabel={ '113°E','114°E', ...
    '115°E', '116°E', '117°E'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:10:92; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=6:30:276; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
[c51,h51] = contourf(distance,wavelength,real_wave,levels,'b--','linewidth',1.2); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('高度:12km');
set(title('高度:12km'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
%colorDepth = 1000;
%colormap(jet(colorDepth));
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
[c52,h52] = contour(distance,wavelength,real_wave,levels,'r-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off

%%%%%%%%%%Figure 6%%%%%%%%%%%%%%%%
subplot(3,2,6);
%%%%%%%%%%%Figure 6%%%%%%%%%%%%%%%
%%%%%%%%%%%%ncread(source,varname)%%%%%%%%%%%%%
distance_lon=ncread('./zll_sunwei_F.hdfw2.nc','lon');
w_velocity=ncread('./zll_sunwei_F.hdfw2.nc','w');
%%%%%%%%%%%ncread(source,varname,start,count,stride)%%%%%%%%%%%%%%%
distance=ncread('./zll_sunwei_F.hdfw2.nc','lon', 65, 136, 1);
w=ncread('./zll_sunwei_F.hdfw2.nc','w',[65,126,7],[136,1,1],[1,1,1]);
%count=[136,1,1]从每一维的start开始读取的总数目%%stride=[1,1,1]设置读取的步长%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n=length(distance);
n=length(w);
%distance =(0:3:3*(n-1)); 
%2%%%%%%%%%%begin variance%%%%%%%%%%%%%%%%%%%%%
%variance = std(w(:))^2;
variance = std(w)^2;
w=(w - mean(w))/sqrt(variance);
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
pad       = 1;           % if set to 1 (default is 0)
dj          = 1/12;      %the spacing between discrete scales. Default is 0.25%
s0          = dstride;  %S0 = the smallest scale of the wavelet.  Default is 2*DT%%
j1          = 6.5/dj;    %Default is J1 = (LOG2(N DT/S0))/DJ%
lag1      = 0.72;        
[wave,wavelength,scale,coi] = wavelet(w(:),dstride,pad,dj,s0,j1,mother);
real_wave = real(wave); 
                            %WAVE is the WAVELET transform of Y. This is a complex array
                            %of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
                            %ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
                            %The WAVELET power spectrum is ABS(WAVE)^2.
                            %Its units are sigma^2 (the time series variance)
%3%%%%%end wavelet transformation%%%%%%%%%%%%%%%%%%%%%%

%4%%%begin Figure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(gcf,'Position', [400,100,300,300]);
%设置绘图的大小,不需要到word里再调整大小(7cm)%
%set(gca,'Position',[.13 .17 .80 .74]); 
%get hanlde to current axis返回当前图形的当前坐标轴的句柄%
%设置xy轴在图片中占的比例，可能需要自己微调%%
%%%%%Figure resize%%%%%%%%%%%
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
%将字体大小改为8号字，在小图里很清晰%%%
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%将线宽改为2%%'Fontname', 'Times newman'%%

%%%%%%%%%%PLOT%%%%%%%%%%%%
levels = (0.025:0.25:5.5); 
%levels = (-5.5:0.065:5.5);  %for color figure 
%v = [0,0.5,1.0,1.5];
Xtick_pos=65:33:201; %确定Label显示的位置
XtickLabel={ '113°E','114°E', ...
    '115°E', '116°E', '117°E'};
%Yticks = fix(min(wavelength)):(max(wavelength)+min(wavelength))/2:fix(max(wavelength));
%Yticks=0:100:490;
Ytick_pos=2:10:92; %确定Label显示的位置;Yticks=0:5:50
YtickLabel=6:30:276; %and modify ylim=[0,60];
%Yticks = fix(min(wavelength)):5:fix(max(wavelength));
[c1,h1] = contourf(distance,wavelength,real_wave,levels,'b--','linewidth',1.2); 
%contour(X,Y,Z,V) draw a contour line for each level specified in vector V.
% contour(X,Y,Z,[v v]) to draw contours for the single level v.
%set(h,'linestyle','none');
%clabel(c1,h1);
title('高度:12km');
set(title('高度:12km'),'FontName','Times New Roman','FontSize',18,'Color','b')
xlabel(' ')
ylabel('wavelength(km)')
set(ylabel('wavelength(km)'),'FontSize',18,'Color','b')
%title(strcat(num2str(height(3)),'hPa'))
xlim = [distance(1),distance(length(distance))];   
%ylim = [2,fix(max(wavelength))];
ylim=[2,92];%设置坐标轴上下限
set(gca,'XLim',xlim(:), ...
    'XTick',Xtick_pos(:), ...
    'XTickLabel', XtickLabel);
%get(gca,'xlim');是获取最大最小刻度的%
set(gca,'YLim',ylim(:), ...
    'YDir','default', ...
'YTick',Ytick_pos(:), ...
'YTickLabel',YtickLabel)
%colorDepth = 1000;
%colormap(jet(colorDepth));
set(gca,'FontName','Times New Roman','FontSize',18);
%坐标轴字体％
%set(gca,'Color',[1,0.4,0.6])
 hold on
 %hcb = colorbar('location','EastOutside');
%levels = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0];
levels = (-5.5:0.25:-0.025);
%v = [-0.5,-1.0,-1.5];
[c2,h2] = contour(distance,wavelength,real_wave,levels,'r-','linewidth',1.2);
%clabel(c2,h2);
%colormap(gca,'jet')
set(gca, 'XColor','blue');       % X轴的颜色%
set(gca, 'YColor','blue');        % Y轴的颜色%
hold off

%%%%%%%%%%%%%%end Figure %%%%%%%%%%%%%%%%%
%5%%%%%%%begin output%%%%%%%%%%%%%%%%%
 %t=strcat('wavelet_w_255',height_12km);
 %print ( t ,'-djpeg ','-r300')；
 %set(gcf, 'color','white','PaperPositionMode', 'auto')   % Use screen size%
 set(gcf, 'color','white','PaperPositionMode', 'manual')  
 set(gcf,'color','white'); %设定figure的背景颜色%
 %set(gcf,'paperunits','inches');
 %set(gcf,'PaperPosition',[100 100 800 600]);
 set(gcf,'unit','centimeters','position',[1 1 50 50]);
 %%%%%将生成的图形保存为7cm*5cm%%%%%%%%%
 %print ( gcf ,'-djpeg ','-r300','wavelet-w.jpg');
 exportfig(gcf,'pic3.eps', 'width',12, 'fontmode','fixed', 'fontsize',18, 'color', 'cmyk');
 %exportfig(gcf, 'fig2.eps', 'FontMode', 'fixed', 'FontSize', 10, 'color', 'cmyk' );
 %t=strcat('wavelet-w');
 %print ( t ,'-djpeg ','-r300');
 %print ( t ,'-dbitmap ','-r300');
 
 %%%指定存储格式。常用的有：%%
 %%%PNG格式：，‘-dpng’（推荐这一种，与bmp格式一样清晰，文件也不大）%%
 %%%JPEG：    ‘-djpeg’（文件小，较清晰）%%
 %%%TIFF：     ‘-dtiff’%%
 %%% BMP：     ‘-dbitmap’（清晰，文件极大）%%
 %%%GIF：      ‘-dgif’   （文件小但不清晰）%%
 %print(gcf,'-dpdf','wavelet-w.pdf')
 %print (gcf,'-dpdf ','wavelet-w.pdf' );
 % 把上面画的图（句柄为 fig ）保存为 fig2.eps, 字号为 10，彩色%%
%5%%%%%%%%%end%%%%%%%%%%%%%%%%%%%%
 