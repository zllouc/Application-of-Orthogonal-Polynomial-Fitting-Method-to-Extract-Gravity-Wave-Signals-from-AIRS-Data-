clear all;clear;clc;
%%%%%%%%%%%%get 4.3?m radiance and TBB
import matlab.io.hdfeos.*
FILE_NAME1='AIRS.2013.07.12.172.L1B.SUBIBRAD.v5.0.22.0.G13311152621.hdf';
SWATH_NAME='L1B_AIRS_Science';
file_id1 = sw.open(FILE_NAME1, 'rdonly');
swath_id1 = sw.attach(file_id1, SWATH_NAME);
% Read data from a data field.
DATAFIELD_NAME='radiances';
long1 = sw.readField(swath_id1, 'Longitude', [], [], []);
lati = sw.readField(swath_id1, 'Latitude', [], [], []);
t = sw.readField(swath_id1, 'Time', [], [], []);
sw.detach(swath_id1);
sw.close(file_id1);

%%%%%%%%%%read 4.3um TBB disturblance
% tb_c=csvread('/Users/nyp/Desktop/nyp_pro/get_data/tb4_3.csv');
% t_cress=csvread('/Users/nyp/Desktop/graduate/smallpaper/polyminal/savedat/tbb_cresssman3_100.csv')
tb_4_3=xlsread('tb_c_4.xls');
t_cress=csvread('tbb_cresssman3_100.csv');
dis_4_3_cress=tb_4_3-t_cress;
t_opf=csvread('tbb_opf67.csv');
dis_4_3_opf=tb_4_3-t_opf;
dis_4_3_5=xlsread('dis1_5.xls');
dis_4_3_4=xlsread('dis1_4.xls');
m_4_3_cress=tb_4_3-dis_4_3_cress;%m_4_3_cress=tb_c-dis_cress;% m_cress=tb_c-dis_cress;
m_4_3_opf=tb_4_3-dis_4_3_opf;
m_4_3_4=tb_4_3-dis_4_3_4;
m_4_3_5=tb_4_3-dis_4_3_5;
% bkbt={m_4_3_4;m_4_3_cress;m_4_3_5;m_4_3_opf};
disbt={dis_4_3_4;dis_4_3_cress;dis_4_3_5;dis_4_3_opf};

% dis_cress=tb_c-t_cress;
for ii=1:4
%%%%%To find the data nearest 21?N
for u=1:90
     for v=1:135
        lon_real=long1;
        lat_real=lati;
        po=find(lat_real>21 & lat_real<22);
        m_21bk{ii}=disbt{ii}(po);
        lat_21=lati(po);
        lon_21=long1(po);
        end
end   
[long_21,id] =sort(lon_21);% sort(x) means from small numer to big
m4_21bk{ii}=m_21bk{ii}(id);%%Arrange the corresponding data according to longitude from large to small
%m_zuo{ii}=m4_21bk{ii}(1:562);
%m_zuo{ii}=m4_21bk{ii}(1:69);%%0-198.9km
%m_you{ii}=m4_21bk{ii}(494:562);%%868.53-1070.9km
end
lati_21=lat_21(id);%%按照经度从大到小将相应的数据排列
%%%%%%%%%%dis4_21 is the TBB disturblance nearest the 21?N 
%%%%%%%%%%Get the approximate number of km for each longitude 
maxlon=max(long_21);minlon=min(long_21);
km=-111*cos(21)*(maxlon-minlon); %1?longitude~111*cos(?)km ?-->latitude 把每个经度大约对应的km数得出
dkm1(1:561)=0;
dkm1=dkm1';
for ii1=2:id
        dkm(ii1)=((long_21(ii1)-long_21(ii1-1))/(maxlon-minlon))*km;
        dkm1(ii1)=dkm1(ii1-1)+dkm(ii1);
end 
tt=['(a)';'(b)';'(c)';'(d)'];
dkmzuo=dkm1(1:562);
%dkmyouorig=dkm1(494:562);
%dkmyou=dkmyouorig-dkmyouorig(1)+dkmzuo(69);
for ij=1:4
%%%%%%%zuo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%left
datcressman=m4_21bk{ii};
variance1 = std(datcressman)^2;
dat_cress = (datcressman - mean(datcressman))/sqrt(variance1) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算采样率为DT的矢量Y（长度为N）的小波变换
% dof=1;
n = length(dat_cress);%数据的矢量长度 The vector length of the data
dt = 1.7184;% 每个Y值之间的时间量
s0 = 2*dt; 
time=dkmzuo*dt;% time = [0:length(data)-1]*dt + 0;  % construct time array
% xlim = [0,dkmzuo(end)];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.001;    % this will do 10 sub-octaves per octave . 
j1 =fix((log(n*dt/s0)/log(2))/dj);% this says start at a scale of 6 months
%% j1 = 10/(dj);    % this says do 7 powers-of-two with dj sub-octaves each
% X=data(1:n-1);Y=data(2:n);
% lag1=corr(X,Y);%%autocorrelation for red noise background
mother = 'Morlet';
lag1=0.72;

save dat_cress -ascii dat_cress
% Wavelet transform:
[wave,period,scale,coi] = wavelet(dat_cress,dt,pad,dj,s0,j1,mother);% scale = s0*2.^((0:J1)*dj);
power = (abs(wave)).^2 ;  %求小波的功率谱(模的平方)  
realpart =real(wave);   %求小波的实部
modulus =abs(wave);   %求小波的模
phase =atan2(imag(wave),real(wave));   %求小波的阶，波位相
variance =sum(power')/n;  %计算小波方差
%%%scale 尺度指数矢量,由S0*2^(j*DJ), j=0...J1给出，在J1+1处是尺度的总
%%%period傅里叶周期（在时间单位）的矢量，与尺度相关。
%%%COI = 如果指定了，就返回 Cone-of-Influence，这是一个N点的矢量，包含了在那特定时间有用信息的最大周期。大于这的周期wie边缘效应。这可以用来在等值线画图上绘制COI线
% Significance levels: (variance=1 for the normalized SST)显著性水平
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:小波全谱&显著性水平
global_ws = variance1*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(variance1,dt,scale,1,lag1,-1,dof,mother);

% Scale-average between El Nino periods of 2--x years
avg = find((scale >= 2) & (scale < 7));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = power ./ scale_avg;   % [Eqn(24)]
scale_avg = variance1*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(variance1,dt,scale,2,lag1,-1,[2,6.9],mother);
%%%%%%%right%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%datyou=m_you{ij};
%variancey = std(datyou)^2;
%dat_you = (datyou - mean(datyou))/sqrt(variancey) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算采样率为DT的矢量Y（长度为N）的小波变换
% dof=1;
%ny = length(dat_you);%数据的矢量长度 The vector length of the data
%dty = 1;% 每个Y值之间的时间量
%s0y = 2*dty; 
%timeyou=dkmyou*dty;% time = [0:length(data)-1]*dt + 0;  % construct time array
% xlim = [0,dkmzuo(end)];  % plotting range
%pady = 1;      % pad the time series with zeroes (recommended)
%djy = 0.001;    % this will do 10 sub-octaves per octave . 
%j1y =fix((log(ny*dty/s0y)/log(2))/djy);% this says start at a scale of 6 months
% j1 = 10/(dj);    % this says do 7 powers-of-two with dj sub-octaves each
% X=data(1:n-1);Y=data(2:n);
% lag1=corr(X,Y);%%autocorrelation for red noise background
%mother = 'Morlet';
%lag1y=0.72;
% Wavelet transform:
%[waveyou,periodyou,scaleyou,coiyou] = wavelet(dat_you,dty,pady,djy,s0y,j1y,mother);% scale = s0*2.^((0:J1)*dj);
%poweryou = (abs(waveyou)).^2 ;  %求小波的功率谱(模的平方)  
%realpartyou =real(waveyou);   %求小波的实部
%modulusyou =abs(waveyou);   %求小波的模
%phaseyou =atan2(imag(waveyou),real(waveyou));   %求小波的阶，波位相
% varianceyou =sum(poweryou')/ny;  %计算小波方差
%%%scale 尺度指数矢量,由S0*2^(j*DJ), j=0...J1给出，在J1+1处是尺度的总
%%%period傅里叶周期（在时间单位）的矢量，与尺度相关。
%%%COI = 如果指定了，就返回 Cone-of-Influence，这是一个N点的矢量，包含了在那特定时间有用信息的最大周期。大于这的周期wie边缘效应。这可以用来在等值线画图上绘制COI线
% Significance levels: (variance=1 for the normalized SST)显著性水平
%[signifyou,fft_theoryou] = wave_signif(1.0,dty,scaleyou,0,lag1y,-1,-1,mother);
%sig95you = (signifyou')*(ones(1,ny));  % expand signif --> (J+1)x(N) array
%sig95you = poweryou ./ sig95you;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:小波全谱&显著性水平
%global_wsyou = variancey*(sum(poweryou')/ny);   % time-average over all times
%dofyou = ny - scaleyou;  % the -scale corrects for padding at edges
%global_signifyou = wave_signif(variancey,dty,scaleyou,1,lag1y,-1,dofyou,mother);

% Scale-average between El Nino periods of 2--x years
%avgyou = find((scaleyou >= 2) & (scaleyou < 7));
%Cdeltayou = 0.776;   % this is for the MORLET wavelet
%scale_avgyou = (scaleyou')*(ones(1,ny));  % expand scale --> (J+1)x(N) array
%scale_avgyou = poweryou ./ scale_avgyou;   % [Eqn(24)]
%scale_avgyou = variancey*djy*dty/Cdeltayou*sum(scale_avgyou(avgyou,:));   % [Eqn(24)]
%scaleavg_signifyou = wave_signif(variancey,dty,scaleyou,2,lag1y,-1,[2,6.9],mother);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

whos
color1=[173/255,174/255,245/255];
%--- Contour plot wavelet power spectrum小波功率谱
figure(410)

subplot(2,2,ij);
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
%contour(time,log2(period),sig95,[-99,1],'r','LineWidth',1);hold on;
%clabel(cy,hy,vv1you,'fontsize',10);hold on;
%%%%%%%%%%%%%%%%%%%
%xlabel('D(km)','fontname','Times New Roman','fontsize',18)
%ylabel('L(km)','fontname','Times New Roman','fontsize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%--- Plot global wavelet spectrum
%subplot('position',[0.77 0.37 0.2 0.28])
%plot(global_ws,log2(period))
plot(log2(period),global_ws,'LineWidth',1)

hold on
%plot(global_signif,log2(period),'--')
plot(log2(period),global_signif,'--','LineWidth',1)
hold off
xlabel('Wavelenth (km)')
ylabel('Spectrum')
%title('Global Wavelet Spectrum')
%set(gca,'YLim',log2([min(period),max(period)]), ...
%	'YDir','reverse', ...
%	'YTick',log2(Yticks(:)), ...
%	'YTickLabel','')
set(gca,'XLim',log2([min(period),max(period)]),'FontSize',18);
%set(gca,'XLim',log2([min(period),max(period)]), ...
%	'XDir','reverse', ...
%	'XTick',log2(Yticks(:)), ...
%	'XTickLabel','')
set(gca,'YLim',[0,1.25*max(global_ws)],'FontSize',18)
%set(gca,'FontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ij==1
title([tt(ij,:),' ','4th PF'],'fontname','Times New Roman','fontsize',18)
end
if ij==2
    title([tt(ij,:),' ','Cresssman'],'fontname','Times New Roman','fontsize',18)
end
if ij==3
    title([tt(ij,:),' ','5th PF'],'fontname','Times New Roman','fontsize',18)
end
if ij==4
    title([tt(ij,:),' ','OPF'],'fontname','Times New Roman','fontsize',18)
end
hold on;

%colorbar
% new_fig_handle = shift_axis_to_origin(figure(410));
end
% set(gcf, 'PaperPositionMode', 'auto')

set(gcf,'color','w');
% print('-dtiff','-r600','E:\假期\savedata\chapter4fig\4_26\4_10wavedis.tiff');
% print('-dmeta','E:\假期\savedata\chapter4fig\4_26\4_10wavedis.emf');
% print 'E:\假期\savedata\chapter4fig\4_26\4_10wavedis.eps' -depsc2 -r600
% saveas(gcf,'E:\假期\savedata\chapter4fig\4_26\4_10wavedis.fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4_11

% figure(2)
% subplot(1,3,2);
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% % contour(time,log2(period),sig95,[-99,1],'r','LineWidth',1);hold on;
% sig=sig95;
% sig(sig95<1 & sig95 >-99)=nan;
% contourf(time,log2(period),sig,'linestyle','none');
% colormap(color1);
% hold on;
% % levels = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9];
% levels = [0.5,1,2,3,4,5,6,7,8,9];
% v = [3,5,7,9];
% % Yticks = 0:50:562;
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% [c,h]=contour(time,log2(period),real(wave),levels,'k-');
% clabel(c,h,v,'fontsize',6);
% xlabel('Distance/km','FontSize',10)
% ylabel('Wavelenghth/km','FontSize',10)
% set(gca,'XLim',xlim(:))
% hold on;
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','default','YTickLabel',Yticks)
% hold on;
% % levels = [-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5,-5,-5.5,-6,-7,-8,-9];
% levels = [-0.5,-1,-2,-3,-4,-5,-6,-7,-8,-9];
% v = [-1,-3,-5,-7,-9];
% [c,h] = contour(time,log2(period),real(wave),levels,'r--');
% clabel(c,h,v,'fontsize',6);
% hold on
% % title('Brightness Temperature Disturbance Wavelet Power Spectrum','FontSize',10)  %图题
% title('(b)','FontSize',10)  %图题
% % cone-of-influence, anything "below" is dubious
% % 锥形影响域(cone of influence，COI),影响锥之下的都是不可信的
% plot(time,log2(coi),'b','LineWidth',1)
% grid on;
% set(gca,'GridLineStyle',':')
% hold off;
% 
% % figure(3)
% % subplot(2,2,1);
% % % plot(log2(period),global_ws)
% % % hold on
% % % plot(log2(period),global_signif,'--')
% % % hold off
% % % ylabel('Power (degC^2)')
% % % Xticks = 0:50:562;
% % % title('(b)')
% % % set(gca,'XLim',log2([min(period),max(period)]), ...
% % % 	'XDir','default', ...
% % % 	'XTick',log2(Xticks(:)), ...
% % % 	'XTickLabel','')
% % % set(gca,'YLim',[0,1.25*max(global_ws)])
% % % hold off
% % plot(global_ws,log2(period),'k')
% % hold on
% % plot(global_signif,log2(period),'r--')
% % hold off
% % xlabel('Power/degC^2')
% % title('(a)','FontSize',10)  %图题
% % set(gca,'YLim',log2([min(period),max(period)]), ...
% % 	'YDir','default','YTickLabel',Yticks)
% % hold on;
% % ylabel('Wavelenghth/km','FontSize',10)
% % hold on;
% % set(gca,'XLim',[0,1.25*max(global_ws)])
% % % %--- Plot global wavelet spectrum小波全谱
% % figure(4);
% % % subplot(2,2,3);
% % plot(log2(period),global_ws)
% % hold on;
% % plot(global_signif,log2(period),'--')
% % hold off;
% % xlabel('Power','FontSize',9.5)
% % title('GWS','FontSize',10)          %图题
% % set(gca,'YLim',log2([min(period),max(period)]), ...
% % 	'YDir','default', ...
% % 	'YTick',log2(Yticks(:)), ...
% % 	'YTickLabel','')
% % set(gca,'ygrid','on');
% % set(gca,'GridLineStyle',':')
% % set(gca,'XLim',[0,1.25*max(global_ws)])
% hold off
% 
% % figure(3)
% % subplot(2,2,1);
% % plot(log2(period),global_ws)
% % hold on
% % plot(log2(period),global_signif,'--')
% % hold off
% % ylabel('Power (degC^2)')
% % Xticks = 0:50:562;
% % Ylimit = 0:50:562;
% % 
% % title('(a)')
% % set(gca,'XLim',xlim(:))
% % hold on;
% % set(gca,'YLim',log2([min(period),max(period)]), ...
% % 	'YDir','default','YTickLabel',Yticks)
% % hold on;
% % 
% % set(gca,'XLim',log2([min(period),max(period)]), ...
% % 	'XDir','default', ...
% % 	'XTick',log2(Xticks(:)), ...
% % 	'XTickLabel','')
% % set(gca,'YLim',[0,1.25*max(global_ws)])
% % 
% % hold off