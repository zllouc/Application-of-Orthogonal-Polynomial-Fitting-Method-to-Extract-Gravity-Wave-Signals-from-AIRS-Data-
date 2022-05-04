clc;clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file2='./wrfout_d02_2013-07-12_12_00_00.nc';
% ncdisp(file2);
times=ncread(file2,'Times');hgt2=ncread(file2,'HGT');
%uwind2=ncread(file2,'U');vwind2=ncread(file2,'V');
%w_d02=ww2(:,:,3,16);ZNW =ncread(file2,'ZNW');
lond2=ncread(file2,'XLONG');latd2=ncread(file2,'XLAT');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='./wrfout_d02_2013-07-12_12_00_00.nc';
[var1,~] = Get_InterpVar( filename,'W',26,16,0,'z' );  
%%%用这个函数插值到指定层次，得到插值后的变量%%%
lon2=lond2(:,:,16);lat2=latd2(:,:,16);
wname=var1;
tit=['(a)'];
hours=['26km'];
%%%%%%%%%%%%%%%%%%%%%%%%
%%skpw2=4;[N,M]=size(lon2);
%%%%%%%%%%%%
latmin = 17;
latmax = 30;
latint = 10;
lonint = 10;
lonmin = 115;
lonmax =132;
yt = [17 23 30];
xt = [115 123 132];
xtl=['115°E';'123°E';'132°E'];
ytl=['17°N';'23°N';'30°N'];
left = [.06,.545];
bottom = [.59,.17];
width = .415;
height = .325;
%%%%%%%%%%
figure(1)
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;
m_coast('color',[0.7 0.7 0.7]);hold on;
dlevels = [-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4] ;
wnamez=makecolor(wname,dlevels);
m_contourf(lon2,lat2,wnamez,'linestyle','none');
hold on;% m_contourf(long1,lati,real(tbb));
m_coast;
hold on;
%colorbar;% 'FontWeight', 'bold',
%colorbar('FontSize',18,'Color','b');
cbwt = .49;
h1= colorbar('FontSize',12);% 'FontWeight', 'bold',
%h1= colorbar('SouthOutside','Position',[.5-cbwt/2,.085,cbwt,.02],'FontSize',12);% 'FontWeight', 'bold',
set(h1,'Ticks',[0,1,2,3,4,5,6,7,8],'TickLabels',dlevels) ;
set(get(h1,'ylabel'),'string','m.s^-^1','fontname','Times New Roman',...
    'fontsize',18);%,[a b]为绘图区域左下点的坐标。c，d分别为绘图区域的宽和高。 

%%m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
%%m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',14,...
%%'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
%%hold on;
%%m_coast('color',[0.7 0.7 0.7]);hold on;
%%dlevels2 = [0,5,10,15,20,25,30,35,40] ;
 
set(gcf, 'PaperPositionMode', 'auto')
print('-dtiff','-r600','./7_2n3.tiff');
saveas(gcf,'./7_2n3.fig');
% print 'E:\假期\savedata\chapter4fig\4_26\6_3tbb2.eps' -depsc2 -r600

% % subplot(3,3,2)
% % 
% % % contourf(w21);  
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % caxis([-0.5,0.5]);
% % hold on;
% % m_contourf(lon12,lat12,w15,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % m_coast;
% % hold on;
% % colormap(othercolor('BuDRd_18'));
% % % colorbar;
% % units = 'K';
% % hold on;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % title('(b)','FontSize',10)
% % 
% % hold on;
% % % figure(5)
% % subplot(3,3,3)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w16,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % m_coast;
% % hold on;
% % colormap(othercolor('BuDRd_18'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(c)','FontSize',10)
% % 
% % hold on;
% % subplot(3,3,4)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w17,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % 
% % % colormap(othercolor('BuGn6'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(d)','FontSize',10)
% % hold on;m_coast;
% % hold on;
% % subplot(3,3,5)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w18,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % 
% % % colormap(othercolor('BuGn6'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(e)','FontSize',10)
% % hold on;m_coast;
% % hold on;
% % subplot(3,3,6)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w19,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % 
% % % colormap(othercolor('BuGn6'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(f)','FontSize',10)
% % hold on;m_coast;
% % hold on;
% % subplot(3,3,7)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w21,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % 
% % % colormap(othercolor('BuGn6'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(g)','FontSize',10)
% % hold on;m_coast;
% % hold on;
% % subplot(3,3,8)
% % m_proj('miller','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
% %        'xaxislocation','bottom','xlabeldir','middle','linestyle','none','ticklen',.02);
% % hold on;
% % caxis([-0.5,0.5]);
% % m_contourf(lon12,lat12,w26,'linestyle','none');
% % hold on;% m_contourf(long1,lati,real(tbb));
% % 
% % % colormap(othercolor('BuGn6'));
% % % colorbar;
% % m_text(122.71,24.98,'o','Color','k','FontSize',9);
% % hold on;
% % m_plot(tclonlat1212_1312(:,2),tclonlat1212_1312(:,1),'LineWidth',1.5,'color','k');
% % hold on;
% % units = 'K';title('(h)','FontSize',10)
% % hold on;m_coast;
% % hold on;
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%剖线


% subplot(1,2,1)
%m_proj('mercator','lat',[17 31],'lon',[115 132]);m_grid('tickdir','out','yaxislocation','left',...
%       'xaxislocation','bottom','xlabeldir','middle','linestyle',':','ticklen',.02);
%hold on;
%m_plot(uslat,uslon,'LineWidth',1,'color','k');%,'color','k'
%hold on;
%m_coast('patch',[.7 .7 .7]);%m_coast;
%hold on;
% title('(a)','FontSize',10)% figure(7)
% subplot(1,2,2)
% tii=1:62;
% plot(tii,usp,'LineWidth',1.5,'color','k');
