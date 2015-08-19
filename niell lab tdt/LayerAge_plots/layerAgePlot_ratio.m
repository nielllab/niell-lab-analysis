function [meandata errdata N]=layerAgePlot_ratio(data1,data2,ageList,layer,inh,used,label,titlestr);
ageList=ageList';
colors = 'bmrg'
    for age = 1:2
        
      for group=1:4
          if group==1
            uselist = (ageList==age & (layer<=3) & ~inh & used');
            elseif group==2
            uselist =(ageList==age & (layer==4) & ~inh & used');
            elseif group==3
            uselist =(ageList==age & (layer==5) & ~inh & used');
            elseif group==4
            uselist =(ageList==age & (layer==6) & ~inh & used');      
           end   
                  
                     if sum(uselist)>2
                         
                     mx=nanmedian(data1(uselist));
                     sx= semedian(data1(uselist));
                     
                     my=nanmedian(data2(uselist));
                     sy=semedian(data2(uselist));
%                      
                     mediandata_x(group,age)=mx;
                     errdata_med_x(group,age)=sx;
                     
                     mediandata_y(group,age)=my;
                     errdata_med_y(group,age)=sy;
                     
                     ratio_ny_nx(group,age)=nanmedian(data2(uselist))/nanmedian(data1(uselist));
                     sem_ratio=semedian_ratio(data2(uselist),data1(uselist))
                     err_ratio(group,age)=sem_ratio
                     
                     
                     end
                     
                    
      end
      end
                  
                    
%                     figure
%                     errorbar(1:4,mediandata_x, errdata_med_x,'k');hold on
%                     errorbar(1:4,mediandata_y, errdata_med_y,'g');hold on
%                     
%                     figure
%                     errorbar(1:4,ratio_ny_nx,err_ratio,'k');hold on
                    
                    figure
                    barweb(ratio_ny_nx,err_ratio)
                    ylabel(label);
                    set(gca,'Xtick',1:5);
                    
       end
             
             
      
    


      

   




    
