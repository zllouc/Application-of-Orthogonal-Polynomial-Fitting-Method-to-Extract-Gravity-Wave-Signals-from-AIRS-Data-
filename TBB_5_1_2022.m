%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
datadir='./';
filelist=dir([datadir,'*.nc']);
mmw=10;nnw=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  set(gcf,'position',[0,0,750,750]);
  %%%%gcf是当前figure，可以改变当前图框的大小
  s1=1;
  filename=[datadir,filelist(s1).name];
  ncid=netcdf.open(filename,'NC_NOWRITE');
  va4  = ncread(filename,'bt_4mu_var'); %读入变量precip
  pt4 = ncread(filename,'bt_4mu_pt');
  bt4=ncread(filename,'bt_4mu');
  bt81  = ncread(filename,'bt_8mu'); 
  va15  = ncread(filename,'bt_15mu_low_var'); %读入变量precip
  pt15 = ncread(filename,'bt_15mu_low_pt');
  bt15=ncread(filename,'bt_15mu_high');
  TimeData  = ncread(filename,'time'); %读入变量time
  lonreal1  = ncread(filename,'lon'); %读入变量lon
  latreal1  = ncread(filename,'lat'); %读入变量lat
  [m,n]=size(lonreal1);
  netcdf.close(ncid); 
  gv4=va4;% 关闭文件%%% nn=150;n2=200; 15:00-20:00 %%nn=20;n2=90; 2:00-9:00
  nn=150;n2=200;nn1=135*nn+1;nn2=135*n2;
  gvbt8_7_9=bt81(1:m,nn1:nn2);
  lon1=lonreal1(1:m,nn1:nn2);lat1=latreal1(1:m,nn1:nn2);
% save('E:\假期\savedata\chapter5\5_1atbb','gvbt8_7_9')
% save('E:\假期\savedata\chapter5\5_1alon','lon1')
% save('E:\假期\savedata\chapter5\5_1alat','lat1')
   clear nn n2 nn1 nn2
  %%%%%%%%%%%%%%%%%%%%%%%%%%7/12/5:00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s2=4;
  filename2=[datadir,filelist(s2).name];
  ncid=netcdf.open(filename2,'NC_NOWRITE');
  bt82  = ncread(filename2,'bt_8mu'); 
  lonreal2  = ncread(filename2,'lon'); %读入变量lon
  latreal2  = ncread(filename2,'lat'); %读入变量lat
  [m,n]=size(lonreal2);
  netcdf.close(ncid); 
  nn=20;n2=90;nn1=135*nn+50;nn2=135*n2+49;
  gvbt8_7_121=bt82(1:m,nn1:nn2);
  lon12=lonreal2(1:m,nn1:nn2);lat12=latreal2(1:m,nn1:nn2);
  clear nn n2 nn1 nn2
% save('E:\假期\savedata\chapter5\5_1btbb.mat','gvbt8_7_121')
% save('E:\假期\savedata\chapter5\5_1blon.mat','lon12')
% save('E:\假期\savedata\chapter5\5_1blat.mat','lat12')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%7/12/17:00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s3=4;
  filename=[datadir,filelist(s3).name];
  ncid=netcdf.open(filename,'NC_NOWRITE');
  bt83  = ncread(filename,'bt_8mu'); 
  lonreal3  = ncread(filename,'lon'); %读入变量lon
  latreal3  = ncread(filename,'lat'); %读入变量lat
  [m,n]=size(lonreal3);
  netcdf.close(ncid); 
  % 关闭文件nn=166;n2=167;
  nn=150;n2=200;nn1=135*nn+1;nn2=135*n2;
  gvbt8_7_122=bt83(1:m,nn1:nn2);
  lon13=lonreal3(1:m,nn1:nn2);lat13=latreal3(1:m,nn1:nn2);
   clear nn n2 nn1 nn2
% save('E:\假期\savedata\chapter5\5_1ctbb.mat','gvbt8_7_122')
% save('E:\假期\savedata\chapter5\5_1clon.mat','lon13')
% save('E:\假期\savedata\chapter5\5_1clat.mat','lat13')
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%7/13/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  s4=5;
  filename=[datadir,filelist(s4).name];
  ncid=netcdf.open(filename,'NC_NOWRITE');
  bt84  = ncread(filename,'bt_8mu'); 
  lonreal4  = ncread(filename,'lon'); %读入变量lon
  latreal4  = ncread(filename,'lat'); %读入变量lat
  [m,n]=size(lonreal4);
  netcdf.close(ncid); 
% 关闭文件
  nn=150;n2=200;nn1=135*nn+1;nn2=135*n2;
  gvbt8_7_13=bt84(1:m,nn1:nn2);
  lon14=lonreal4(1:m,nn1:nn2);lat14=latreal4(1:m,nn1:nn2);
% save('E:\假期\savedata\chapter5\5_1dtbb.mat','gvbt8_7_13')
% save('E:\假期\savedata\chapter5\5_1dlon.mat','lon14')
% save('E:\假期\savedata\chapter5\5_1dlat.mat','lat14')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latmin = 0;
latmax = 40;
latint = 10;
lonint = 10;
lonmin = 100;
lonmax =150;
xt = [lonmin:lonint:lonmax];
yt = [latmin:latint:latmax];
% for i = 1:length(xt)
%     bl = num2str(xt(i));
%    xtl(i,:) = [bl,'°E'];
%    
% end
ytl=['0°  ';'10°N';'20°N';'30°N';'40°N'];
xt=[100 120 140];xtl=['100°E';'120°E';'140°E'];

%%
left = [.06,.545];
bottom = [.57,.11];
width = .415;
height = .345;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
subplot(2,2,1);
i1=1;
subplot('Position',[left(1+mod(i1-1,2)),bottom(ceil(i1/2))+0.05,width,height]);
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;
caxis([190,290]);
m_contourf(lon1,lat1,gvbt8_7_9,'linestyle','none');
hold on;% m_contourf(long1,lati,real(tbb));
m_coast('color',[0.7 0.7 0.7]);
hold on;
colormap(othercolor('BuDRd_18'));
aa=get(gca);
xx=aa.XLim;%获取横坐标上下限
yy=aa.YLim;%获取纵坐标上下限
kk1=[0.67 0.18];kk2=[0.67 0.075];%给定text相对位置
x0=xx(1)+kk1(1)*(xx(2)-xx(1));%获取text横坐标
y0=yy(1)+kk1(2)*(yy(2)-yy(1));%获取text纵坐标text(x0,y0,'tex');
x1=xx(1)+kk2(1)*(xx(2)-xx(1));%获取text横坐标
y1=yy(1)+kk2(2)*(yy(2)-yy(1));%获取text纵坐标
ss=[max(max(gvbt8_7_9)),min(min(gvbt8_7_9))];
%text(x0,y0,['T_m_a_x = ',num2str(ss(1),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
%text(x1,y1,['T_m_i_n = ',num2str(ss(2),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
title('(a) 16:42UTC July 9','fontname','Times New Roman','fontsize',18)
hold on;clear ss ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2);
i2=2;
subplot('Position',[left(1+mod(i2-1,2)),bottom(ceil(i2/2))+0.05,width,height]);
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;
caxis([190,290]);
m_contourf(lon12,lat12,gvbt8_7_121,'linestyle','none');
hold on;% m_contourf(long1,lati,real(tbb));
m_coast('color',[0.7 0.7 0.7]);
hold on;
colormap(othercolor('BuDRd_18'));
ss=[max(max(gvbt8_7_121)),min(min(gvbt8_7_121))];
%text(x0,y0,['T_m_a_x = ',num2str(ss(1),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
%text(x1,y1,['T_m_i_n = ',num2str(ss(2),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
title('(b) 05:00UTC July 12','fontname','Times New Roman','fontsize',18)
hold on;clear ss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
i3=3;
subplot('Position',[left(1+mod(i3-1,2)),bottom(ceil(i3/2))+0.05,width,height]);
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;
caxis([190,290]);
m_contourf(lon13,lat13,gvbt8_7_122,'linestyle','none');
hold on;% m_contourf(long1,lati,real(tbb));
m_coast('color',[0.7 0.7 0.7]);
hold on;
colormap(othercolor('BuDRd_18'));
ss=[max(max(gvbt8_7_122)),min(min(gvbt8_7_122))];
%text(x0,y0,['T_m_a_x = ',num2str(ss(1),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
%text(x1,y1,['T_m_i_n = ',num2str(ss(2),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
title('(c) 17:12UTC July 12','fontname','Times New Roman','fontsize',18)
hold on;clear ss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4);
i4=4;
subplot('Position',[left(1+mod(i4-1,2)),bottom(ceil(i4/2))+0.05,width,height]);
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;
caxis([190,290]);
m_contourf(lon14,lat14,gvbt8_7_13,'linestyle','none');
hold on;% m_contourf(long1,lati,real(tbb));
m_coast('color',[0.7 0.7 0.7]);
hold on;
colormap(othercolor('BuDRd_18'));
ss=[max(max(gvbt8_7_13)),min(min(gvbt8_7_13))];
%text(x0,y0,['T_m_a_x = ',num2str(ss(1),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
%text(x1,y1,['T_m_i_n = ',num2str(ss(2),'%.3f')],'FontSize',9,'BackgroundColor',[1 1 1]);
title('(d) 17:54UTC July 13','fontname','Times New Roman','fontsize',18)
hold on;clear ss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555
P=get(gca,'Position');cbwt = .49;
h1= colorbar('SouthOutside','Position',[.5-cbwt/2,.085,cbwt,.02],'FontSize',18);% 'FontWeight', 'bold',
set(get(h1,'ylabel'),'string','BT(K)','fontname','Times New Roman',...
    'fontsize',18);%,[a b]为绘图区域左下点的坐标。c，d分别为绘图区域的宽和高。 
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'color','w');
% print 'E:\假期\savedata\chapter4fig\4_26\4_2tbb2.emf' %-dmeta -r300
%print('-dtiff','-r600','5_1tbb2.tiff');
%print('-dmeta','5_1tbb2.emf');
% print 'E:\假期\savedata\chapter4fig\4_26\4_3tbb2.eps' -depsc2 -r600
%saveas(gcf,'5_1tbb2.fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(gcf,'zll-convetion-TBB-5-1.eps','-depsc2', '-r600');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%end of code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

