function va4=get_va(tb,long1,lati)
dr=1;
for u=1:90
     for v=1:135

        lon_real=long1(u,v);
        lat_real=lati(u,v);
        po=find(long1>lon_real-dr & long1<lon_real+dr & lati>lat_real-dr & lati< lat_real+dr);
        tb_part=tb(po);
        tb_ave(u,v)=sum(sum(tb_part));
        [n,c]=size(tb_part);
        nn(u,v)=n*c;

        clear tb_part po
        end
end   
 tb_ave=tb_ave./nn;
tb_gi=( tb-tb_ave ).^2;

  for u=1:90
     for v=1:135
         lon_real=long1(u,v);
         lat_real=lati(u,v);
         po=find(long1>lon_real-dr & long1<lon_real+dr & lati>lat_real-dr & lati< lat_real+dr);
         va4(u,v)=sum(sum(tb_gi(po)));     
     end
  end
va4=va4./(nn-1);
end