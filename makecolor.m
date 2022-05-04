function [z]=makecolor(dismean,dlevels)
z=dismean;
% dlevels = [0,50,200,500,800,1000,1500,2000,3000,4000] ;
% dlevels = [-0.6,-0.5,-0.4,-0.2,0,0.2,0.4,0.5,0.6] ;
 ncolor1=[20 43 140
    0 114 189
    129 184 224
    222 235 250
    243 222 187
    248 118 100
    255 0 0
    204 0 0]/255;
  
for k = 1 : length(dlevels) - 1
     
   z(find(dismean>dlevels(k) & dismean<=dlevels(k+1))) = k-1 ;
   
end
   z(find(dismean<=dlevels(1))) = 0 ;
      z(find(dismean>dlevels(length(dlevels) ))) = length(dlevels)-1 ;

%    cmap = colormap(mycolor2(length(dlevels) - 1)) ;
   cmap = colormap(ncolor1) ;  
   colormap(cmap) ;
    
   caxis([0 length(dlevels)-1]) ;
    
%    cbar = colorbar ;
    
%    set(cbar,'Ticks',[0,1,2,3,4,5,6,7,8],'TickLabels',dlevels) ;
      
hold on
end