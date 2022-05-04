%  GPIdata=GPI_m_w10;  tcdata=tcmesg7_11; color=colorbar;  % 用法1，需注释4行
%  serie_corrcoef(GPI_m_w10,tcmesg7_11,colorbar);          % 用法2 运行此行   

 function [posit,R,P]=serie_corrcoef(GPIdata,tcdata,color)
 % 
monthbegin=7;
monthend=11;
%  热带气旋数据
%  pathx=['E:\work\TCcount\RCSM\tc',num2str(monthbegin),'_',num2str(monthend),'.mat'];
%  tcdata=cell2mat(struct2cell(load(pathx,['tcmesg',num2str(monthbegin),'_',num2str(monthend)])));

% ——————区域说明————————
%     整区 5-30N,100-180E
%     南海 5-30N/5-25N,100-120E
%     西太 5-30N ,  120-180E
%    主要生成区（MDR）：5-30N,105-170E/100-170E

wnplo=21:81;  wnpla=6:31;  % 西太 
scslo=1:21;   scsla=6:26;   % 南海
wslo=1:81;   wsla=6:31;   % 西太-南海
%
%--------------------tc时间序列---------------------------------------------

x(:,1)=mean(tcdata(:,monthbegin:monthend),2);

%--------------------空间场------------------------------------------------
[lo,la,ntm]=size(GPIdata);
nmonth=length(monthbegin:monthend);
ydata=zeros(lo,la,36);
if(ntm>36)
    tm=ntm/nmonth;
    for tt=1:tm
      ydata(:,:,tt)=nansum(GPIdata(:,:,(tt-1)*nmonth+1:(tt-1)*nmonth+nmonth),3);  
    end
elseif(ntm==36)
    tm=ntm;
    ydata=GPIdata;
end
%-------------------------------------------------------------------------

%  求相关
R=zeros(lo,la);
P=zeros(lo,la);
point=1; % 平滑的点数
for i=1:lo
    for j=1:la
       x(:,1)=smooth(x(:,1),point);
       ydata(i,j,:)=smooth(ydata(i,j,:),point);
       [r,p]=corrcoef( x(:,1),ydata(i,j,:) );
       R(i,j) = r(2,1);
       P(i,j) = p(2,1);
    end
end
R=flip(R,2);
P=flip(P,2);

RR=R;
PP=P;
% 判断R值是否大于0.34 小于0.34的值赋值为NaN
rpo=find(abs(RR)<=0.34); %95显著性水平
  RR(rpo)=NaN;
 ppo=find(PP>0.05); % 95显著性水平
  PP(ppo)=NaN;
  
 RR=flip(RR,2);
 PP=flip(PP,2);
 
 r_ave=0;
 if(r_ave==1)
 rwnp=nanmean(nanmean(RR(wnplo,wnpla)));
 rscs=nanmean(nanmean(RR(scslo,scsla)));
 rws=nanmean(nanmean(RR(wslo,wsla)));
 disp('西北太平洋区域平均相关系数');
 disp(num2str(rwnp))
 disp('南海区域平均相关系数');
 disp(num2str(rscs))
 disp('南海-西北太平洋区域平均相关系数');
 disp(num2str(rws))
 end
%% 显著性区域打点
t=1;
for ii=2:81
    for jj=2:41
%         if((abs(R(ii,jj))>=0.34) && (ii~=1)&&(ii~=81)&&(jj~=1)&&(jj~=41))
          if(abs(R(ii,jj))>=0.34)
            posit(t,1)=100+ii-1;
            posit(t,2)=jj-1;
%             posit(t,2)=50-jj+1;
        end
        t=t+1;
    end
end

%% 画图
figure
% colorbar 
%  colorpath='E:\work\colorbar\';
%  color=cell2mat(struct2cell(load([colorpath,'redyellowblue.mat'])));

 % 选取部分 
 coe=RR(:,1:41);  % 0-40N
 cop=PP(:,1:41);
 [lat,lon]=meshgrid(0:1:40,100:1:180);
 m_proj('miller','lon',[100 180],'lat',[0 40]);
 m_contourf(lon,lat,R(:,1:41),'linestyle','none','levelstep',0.01);
 colormap(color);
 m_coast('color','k');
 m_grid('xtick',5,'ytick',5,'linestyle','none','fontsize',12);
 
 hold on
 m_plot(posit(:,1),posit(:,2),'k.','markersize',4);

 title('Correlation between OBS and GPI', 'fontsize',14,'fontname','Times New Roman')
 set(gca,'linewi',1,'fontsize',14);



