clc;clear;
import matlab.io.hdfeos.*
FILE_NAME1='AIRS.2013.07.12.172.L1B.SUBIBRAD.v5.0.22.0.G13311152621.hdf';
SWATH_NAME='L1B_AIRS_Science';
file_id1 = sw.open(FILE_NAME1, 'rdonly');
swath_id1 = sw.attach(file_id1, SWATH_NAME);
% Read data from a data field.
% DATAFIELD_NAME='Geolocation Filed';
DATAFIELD_NAME='radiances';
data1 = sw.readField(swath_id1, DATAFIELD_NAME, [], [],[]);
lon1 = sw.readField(swath_id1, 'Longitude', [], [], []);
lat1 = sw.readField(swath_id1, 'Latitude', [], [], []);
time = sw.readField(swath_id1, 'Time', [], [], []);
scanang1= sw.readField(swath_id1, 'scanang', [], [], []);
% Replace the filled value with NaN.
data11(1:26,:,:)=data1(262:287,:,:);%2322.622-2345.93
data11(27:42,:,:)=data1(294:309,:,:);%2352.54-2366.84
data22=data1(75,:,:);%15?m667.782
data55=data1(160,:,:);%8.1?m，1231wavenumber
% Detach from the Swath object and close the file.
sw.detach(swath_id1);
sw.close(file_id1);
datamean1=mean(data11(:,:,:),1);
datamean2=mean(data22(:,:,:),1);
datamean5=mean(data55(:,:,:),1);
ra1 =squeeze(datamean1);%4.3radiance
ra2 =squeeze(datamean2);%15radiance
ra5=squeeze(datamean5);%8.1radiance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TBB,Plunk
c1 = 1.191066e-5;
c2 = 1.438833;
la1 = 4.3;
la1 = 1e4/la1;
la2 = 15;
la2 = 1e4/la2;
la5=8.1;
la5 = 1e4/la5;
tb_c= (c2*la1)./(log(1+(c1*(la1^3)./ra1)));%4.3Tbb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fig4_3
tbb_poly4= csvread('4_3_a.csv');
tbb_cresssman3_100 = csvread('4_3_b.csv');
dis5=xlsread('4_4_c.xls');
 tbb_opf67= csvread('4_3_d.csv');
dis3_100=tb_c-tbb_cresssman3_100;
dis_poly4=tb_c-tbb_poly4;dis_opf67=tb_c-tbb_opf67;
tbb_poly5=tb_c-dis5;
tti1=['(a)';'(b)';'(c)';'(d)'];
latmin = 0;
latmax = 40;
latint = 10;
lonint = 10;
lonmin = 110;
lonmax =150;
xt = [lonmin:lonint:lonmax];
yt = [latmin:latint:latmax];
for i = 1:length(xt)
    bl = num2str(xt(i));
   xtl(i,:) = [bl,'°E'];
   
end
ytl=['0°  ';'10°N';'20°N';'30°N';'40°N'];
%%
left = [.06,.545];
bottom = [.57,.11];
width = .415;
height = .345;
%%
% xlswrite('4_3_c',tbb_poly5);xlswrite('4_3_d',tbb_opf67); xlswrite('4_3_a',tbb_poly4);xlswrite('4_3_b',tbb_cresssman3_100);
% xlswrite('4_4_c',dis_poly5);xlswrite('4_4_d',dis_opf67); xlswrite('4_4_a',dis_poly4);xlswrite('4_4_b',dis3_100);
fig1={tbb_poly4;tbb_cresssman3_100;tbb_poly5;tbb_opf67};

fig2={dis_poly4;dis3_100;dis5;dis_opf67};

dis_opf67new = dis_opf67(45,:);

fid = fopen('dis_opf67new.txt','wt');
fprintf(fid,'%g\n',dis_opf67new);       % \n 换行
fclose(fid);
%save dis_opf67.txt -ascii dis_opf67
%dis_opf = dis_opf67(45,:);
%save('dis_opf67.txt','dis_opf67','-ascii')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fig4_4
figure(44)
%save dis_poly4 -ascii dis_poly4
% set(gcf,'Units','centimeters','Position',[2.8 0.001 19 15]);%%一张图
% set(gca,'Position',[.0425 .08 .85 .85]);

set(gcf,'Units','centimeters','Position',[2.8 0.001 28 22]);%%4张子图
set(gca,'Position',[.0425 .08 .85 .85]);

for i2=1:4

subplot('Position',[left(1+mod(i2-1,2)),bottom(ceil(i2/2))+0.05,width,height]);
m_proj('Equidistant cylindrical','lon',[lonmin,lonmax],'lat',[latmin,latmax]);
m_grid('xtick',xt,'xticklabel',xtl,'ytick',yt,'yticklabel',ytl,'linestyle','none','fontname','Times New Roman','fontsize',18,...
'tickdir','out','linewidth',.5,'ticklen',.02);% for j = 1:length(yt)
hold on;caxis([-1,1]);

% dlevel2 = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8] ;
% fig2n{i2}=makecolor(fig2{i2},dlevel2);('color',[0.7 0.7 0.7])
m_contourf(lon1,lat1,fig2{i2},'linestyle','none');
%colormap(othercolor('BuDRd_18'));
 
hold on;%
m_coast;
title(tti1(i2,:) ,'fontname','Times New Roman','fontsize',20)%,,'position',[0.1,0.95]
end
P=get(gca,'Position');cbwt = .49;
h1= colorbar('SouthOutside','Position',[.5-cbwt/2,.085,cbwt,.02],'FontSize',18);% 'FontWeight', 'bold',
% set(h1,'Ticks',[0,1,2,3,4,5,6,7,8],'TickLabels',dlevel2) ;
set(get(h1,'ylabel'),'string','Brightness temperature perturbations (K)','fontname','Times New Roman',...
    'fontsize',18);%,[a b]为绘图区域左下点的坐标。c，d分别为绘图区域的宽和高。 
% set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'color','w');
saveas(gcf,'4_4disnew.fig');
% print('-dtiff','-r600','E:\假期\savedata\chapter4fig\4_4disnew.tiff');
% % print('-dmeta','E:\假期\savedata\chapter4fig\4_26\4_4dis.emf');
%print '4_4disnew.eps' -depsc2 -r600;

