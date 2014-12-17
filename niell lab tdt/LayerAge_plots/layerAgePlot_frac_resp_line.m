function [P_FF P_HSF prct_err_FF prct_err_HSF data N ]=layerAgePlot_SF_line(data,ageList,layer,inh,used,used1,label);
ageList=ageList';
figure
hold on

colorlist='bgkmrc';


for age =1:4

    for group=1:5
    
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
            uselist1=(ageList==age & (layer<=3) & ~inh & used1)
        elseif group==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
            uselist1=(ageList==age & (layer==4) & ~inh & used1)
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
            uselist1=(ageList==age & (layer==5) & ~inh & used1)
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
            uselist1=(ageList==age & (layer==6) & ~inh & used1)
        elseif group==5
             uselist = (ageList==age & inh & used);
             uselist1=(ageList==age & inh & used1)
        end


if sum(uselist)>2 
   

    
          
          total(group,age) = sum(~isnan(data(uselist)));
          resp(group,age) = sum(~isnan(data(uselist1)));
          
          frac(group,age)= resp(group,age)/total(group,age);
%           P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
          [M,V]= binostat(total(group,age),frac(group,age));
          
          errdata(group,age) =sqrt(V)/sqrt(total(group,age));
          prct_err(group,age)= errdata(group,age)/resp(group,age);
            
          prct_err_lin(group,age)=prct_err(group,age)*frac(group,age);
end
    end    

end

figure

for group = 1:5
 errorbar(1:4,P_FF(group,:),prct_err_FF(group,:),'color',colorlist(group),'LineWidth',2);
    hold on;
end
 errorbar(3:4,P_FF(5,3:4),prct_err_FF(5,3:4),'color','r','LineWidth',2);
 
 ylabel(label{1,1});
 set(gca,'Xtick',1:4);
 set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});
   

 
 
    end
%         clear HSF M V M1 V1 x f
        
 








    
