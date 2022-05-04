clear;clc;
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
%%%%%%%%%%read 4.3?m TBB disturblance
tb_4_3=xlsread('tb_c_4.xls');
t_cress=csvread('tbb_cresssman3_100.csv');
dis_4_3_cress=tb_4_3-t_cress;
t_opf=csvread('tbb_opf67.csv');
dis_4_3_opf=tb_4_3-t_opf;
dis_4_3_5=xlsread('dis1_5.xls');
dis_4_3_4=xlsread('dis1_4.xls');
m_4_3_cress=tb_4_3-dis_4_3_cress;
m_4_3_opf=tb_4_3-dis_4_3_opf;
m_4_3_4=tb_4_3-dis_4_3_4;
m_4_3_5=tb_4_3-dis_4_3_5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4_6
Abserropoly4=abs(m_4_3_4-tb_4_3);Abserrocress=abs(m_4_3_cress-tb_4_3);
Abserropoly5=abs(m_4_3_5-tb_4_3);Abserroopf=abs(m_4_3_opf-tb_4_3);
% x=0:1.8;y= 
tb4=m_4_3_4(:);tbcress=m_4_3_cress(:);tb5=m_4_3_5(:);tbopf=m_4_3_opf(:);
abs4=Abserropoly4(:);abscress=Abserrocress(:);
abs5=Abserropoly5(:);absopf=Abserroopf(:);
% figure(46)
% set(gcf,'Units','centimeters','Position',[2.8 0 28 35]);
% scatter(abs5,tbopf,8,'r','filled'); hold on;
% scatter(abs5,tb5+6.5,8,[0 0.5 1],'filled');hold on;
% scatter(abscress,tbcress+12.5,8,'g','filled');hold on;
% scatter(abs4,tb4+18.5,8,[0 1 1],'filled');hold on;
% ylim([235.6,262]);
% set(gca,'ytick',[238 245.3 251 257.5 ]);
% set(gca,'yticklabel',{'OPF','5th PF','Cresssman','4th PF'},'fontname','Times New Roman','fontsize',14);
% xlabel('Absolute errors/K','fontname','Times New Roman','fontsize',15)
% % ylabel('Methods','fontname','Times New Roman','fontsize',13)
% % set(gcf, 'PaperPositionMode', 'auto')
% % print('-dtiff','-r600','E:\假期\savedata\chapter4fig\4_26\4_6tbb3.tiff');
% % print('-dmeta','E:\假期\savedata\chapter4fig\4_26\4_6tbb3.emf');
% % % print 'E:\假期\savedata\chapter4fig\4_26\4_3tbb2.eps' -depsc2 -r600
% % saveas(gcf,'E:\假期\savedata\chapter4fig\4_26\4_6tbb3.fig');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4_7
% tb=tb_4_3(:);
% figure(47)
% set(gcf,'Units','centimeters','Position',[2.8 0 28 35]);
% subplot(2,2,1)
% H1=scatter(tb,tb4,8,[0 1 1],'filled');hold on;
% xlim([235.5,243]);
% ylim([235.5,243]);
% plot(tb,tb,'color','k');  hold on;
% legend(H1,'4th PF','location','southeast');hold on;
% xlabel('Estimated brightness temperature/K','fontname','Times New Roman','fontsize',14)
% ylabel('Brightness temperature/K','fontname','Times New Roman','fontsize',14)
% hold on;title('(a)','FontSize',14,'fontname','Times New Roman')
% subplot(2,2,2)
% H2=scatter(tb,tbcress,8,'g','filled');hold on;
% xlim([235.5,243]);ylim([235.5,243]);
% plot(tb,tb,'color','k');  hold on;
% legend(H2,'Cresssman','location','southeast')
% xlabel('Estimated brightness temperature/K','fontname','Times New Roman','fontsize',14)
% ylabel('Brightness temperature/K','fontname','Times New Roman','fontsize',14)
% hold on;
% title('(b)','FontSize',14,'fontname','Times New Roman')
% 
% subplot(2,2,3)
% H3=scatter(tb,tb5,8,[0 0.5 1],'filled');hold on;
% xlim([235.5,243]);ylim([235.5,243]);
% plot(tb,tb,'color','k');  hold on;
% legend(H3,'5th PF','location','southeast')
% xlabel('Estimated brightness temperature/K','fontname','Times New Roman','fontsize',14)
% ylabel('Brightness temperature/K','fontname','Times New Roman','fontsize',14)
% title('(c)','FontSize',14,'fontname','Times New Roman')
% 
% subplot(2,2,4)
% H4=scatter(tb,tbopf,8,'r','filled'); hold on;
% xlim([235.5,243]);ylim([235.5,243]);
% plot(tb,tb,'color','k');  hold on;
% legend(H4,'OPF','location','southeast')
% xlabel('Estimated brightness temperature/K','fontname','Times New Roman','fontsize',14)
% ylabel('Brightness temperature/K','fontname','Times New Roman','fontsize',14)
% title('(d)','FontSize',14,'fontname','Times New Roman')
% % set(gcf, 'PaperPositionMode', 'auto')
% % print('-dtiff','-r600','E:\假期\savedata\chapter4fig\4_26\4_7tbb3.tiff');
% % print('-dmeta','E:\假期\savedata\chapter4fig\4_26\4_7tbb3.emf');
% % % print 'E:\假期\savedata\chapter4fig\4_26\4_3tbb2.eps' -depsc2 -r600
% % saveas(gcf,'E:\假期\savedata\chapter4fig\4_26\4_7tbb3.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%To find the data nearest 21?N
for u=1:90
     for v=1:135
        lon_real=long1;
        lat_real=lati;
        po=find(lat_real>21 & lat_real<22);
        tb_21=tb_4_3(po);
        dis_21_opf=dis_4_3_opf(po);
        dis_21_cress=dis_4_3_cress(po);
        dis_21_4=dis_4_3_4(po);
        dis_21_5=dis_4_3_5(po);
        m_21_cress=m_4_3_cress(po);
        m_21_opf=m_4_3_opf(po);
        m_21_4=m_4_3_4(po);
        m_21_5=m_4_3_5(po);
        lat_21=lati(po);
        lon_21=long1(po);
        end
end   
[long_21,id] =sort(lon_21);% sort(x) means from small numer to big
tb4_21=tb_21(id);
dis4_21_cress=dis_21_cress(id);
dis4_21_opf=dis_21_opf(id);
dis4_21_4=dis_21_4(id);
dis4_21_5=dis_21_5(id);
m4_21_cress=m_21_cress(id);
m4_21_opf=m_21_opf(id);
m4_21_4=m_21_4(id);m4_21_5=m_21_5(id);
%%Arrange the corresponding data according to longitude from large to small
lati_21=lat_21(id);%%按照经度从大到小将相应的数据排列
%%%%%%%%%%dis4_21 is the TBB disturblance nearest the 21?N 
%%%%%%%%%%Get the approximate number of km for each longitude 
maxlon=max(long_21);minlon=min(long_21);
km=-111*cos(21)*(maxlon-minlon); %1?longitude~111*cos(?)km ?-->latitude 把每个经度大约对应的km数得出
dkm1(1:561)=0;
dkm1=dkm1';
for ii=2:id
        dkm(ii)=((long_21(ii)-long_21(ii-1))/(maxlon-minlon))*km;
        dkm1(ii)=dkm1(ii-1)+dkm(ii);
end 
%%%%%%%%%%%%dkm1 is the corresponding distance(km) from the first point
% lon1=long1(:);lat1=lati(:);dis=dis_4_3_5(:);
% dis4_3=[lon1,lat1,dis];
data0=tb4_21;
datacress=dis4_21_cress;
dataopf=dis4_21_opf;
datapoly4=dis4_21_4;datapoly5=dis4_21_5;
data4=dis4_21_5;
data_cress=m4_21_cress;
data_opf=m4_21_opf;
data_4=m4_21_4;
data_5=m4_21_5;
% tb4_km=[data0,dkm1];
% bk2_km=[data_2,dkm1];bk3_km=[data_3,dkm1];bk4_km=[data_4,dkm1];bk5_km=[data_5,dkm1];
% dis4_2_km=[data1,dkm1];
% dis4_3_km=[data2,dkm1];dis4_4_km=[data3,dkm1];dis4_5_km=[data4,dkm1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4_8
figure(49)
set(gcf,'Units','centimeters','Position',[2.8 0 28 35]);ms=20;
%%[10 5 7 5]%图形起点坐标为（10cm，5cm）表示左下点离显示器左侧边界10cm，离下侧边界5cm。
subplot(2,1,1);
% contour(time,log2(period),sig95,[-99,1],'r','LineWidth',1);hold on;
plot(dkm1,data0,'k','LineWidth',2);   hold on;
plot(dkm1,data_4,'color',[0 1 1]);  hold on;
plot(dkm1,data_cress,'color','g');  hold on;
plot(dkm1,data_5,'color',[0 0.5 1]);hold on;
plot(dkm1,data_opf,'r','LineWidth',2);  hold on;
% plot(dkm1,data_4,'color',[1 0 1]);  hold on;
xlabel('Distance/km','fontname','Times New Roman','fontsize',20)
ylabel('Brightness temperature/K','fontname','Times New Roman','fontsize',20)
hold on;
% title('Brightness Temperature Disturbance Wavelet Power Spectrum','FontSize',10)  %图题
title('(a)','fontname','Times New Roman','fontsize',20)
 %set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','default','YTickLabel',Yticks)
hold on;
legend('BT','4th PF','Cressman','5th PF','OPF'...
    ,'location','southeast');hold on;
% figure(2)
subplot(2,1,2);
plot(dkm1,datapoly4,'color',[0 1 1],'LineWidth',2);  hold on;
plot(dkm1,datacress,'g');  hold on;
plot(dkm1,datapoly5,'color',[0 0.5 1]);   hold on;
plot(dkm1,dataopf,'r','LineWidth',2);  hold on;
save datapoly5 -ascii datapoly5
save datapoly4 -ascii datapoly4
save dataopf -ascii dataopf
save datacress -ascii datacress

% plot(dkm1,data4,'color',[1 0 1]);  hold on;
xlabel('Distance/km','fontname','Times New Roman','fontsize',20);
ylabel('Brightness temperature perturbations/K','fontname','Times New Roman','fontsize',20); hold on;
title('(b)','fontname','Times New Roman','fontsize',20)
legend('4th PF','Cressman','5th PF','OPF'...
    ,'location','southeast'); hold on;
set(gcf,'color','w');
colormap(turbo)
freezeColors
hold on;
% set(gcf, 'PaperPositionMode', 'auto')
% print('-dtiff','-r600','E:\假期\savedata\chapter4fig\4_26\4_8tbb3.tiff');
% print('-dmeta','E:\假期\savedata\chapter4fig\4_26\4_8tbb3.emf');
% % print 'E:\假期\savedata\chapter4fig\4_26\4_3tbb2.eps' -depsc2 -r600
% saveas(gcf,'E:\假期\savedata\chapter4fig\4_26\4_8tbb3.fig');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%