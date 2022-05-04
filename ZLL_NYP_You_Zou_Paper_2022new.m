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
load 'you_5pf.dat';
%%load 'datacress.dat';
%%wwa=dataopf;
wwa=you_5pf;
%%wwa=datacress;
%%ncinfo(dataname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end w-velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count=[136,1,1]��ÿһά��start��ʼ��ȡ������Ŀ%%stride=[1,1,1]���ö�ȡ�Ĳ���%%%
%1%%%%%%%%%%%%end ncread%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n=length(distance);
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
dj          =  1/12;      %the spacing between discrete scales. Default is 0.25%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v = [0,0.5,1.0,1.5];
Xtick_pos=1:80:nlength; %ȷ��Label��ʾ��λ��
XtickLabel={'115��E','118��E','121��E','124��E','127��E','130��E','133��E'};
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
subplot(1,2,1);
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
set(title('c) WS(Morlet)'),'FontName','Times New Roman','FontSize',12,'Color','b')
colorDepth = 1000;
colormap(jet(colorDepth));
xlabel('Distance (km)')
ylabel('wavelength(km)')
title('c)WPS(Morlet)')
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
subplot(1,2,2);
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
set(title('Global WPS'),'FontName','Times New Roman','FontSize',12,'Color','b')
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
print(gcf,'zll-you5pf-wavelet.eps','-depsc2', '-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end of code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
