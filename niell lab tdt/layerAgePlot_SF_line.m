function [P_FF P_HSF prct_err_FF prct_err_HSF data N ]=layerAgePlot_SF_line(data,ageList,layer,inh,used,label);
ageList=ageList';
figure
hold on

colorlist='bmkg';

figure

for group=1:5
    for age =1:4
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
        elseif group==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==5
             uselist = (ageList==age & inh & used);
        end


if sum(uselist)>2 
   

    [f,x]=hist(data(uselist),0:0.02:.35);%# create histogram from the data
          
          countdata_FF(group,age) = histc(data(uselist),0.005); %%to use with SF_pref
          N(group,age) = sum(~isnan(data(uselist)));
% % % %         
          P_FF(group,age)= countdata_FF(group,age)/N(group,age)
         
          [M,V]= binostat(N(group,age),P_FF(group,age));
          
          
          errdata_FF(group,age) =sqrt(V)/sqrt(N(group,age));
          prct_err(group,age)= errdata_FF(group,age)/countdata_FF(group,age);
          
          prct_err_FF(group,age)=prct_err(group,age)*P_FF(group,age);
% % % %         
          HSF = histc(data(uselist),0.3:0.05:0.35);
          countdata_HSF(group,age)= HSF(1,1);
% % % %         
          P_HSF(group,age)= countdata_HSF(group,age)/N(group,age)
       
          [M1,V1]= binostat(N(group,age),P_HSF(group,age));
          
          errdata_HSF(group,age) =sqrt(V1)/sqrt(N(group,age));
          prct_err_H(group,age)= errdata_HSF(group,age)/countdata_HSF(group,age);
          
          prct_err_HSF(group,age)=prct_err_H(group,age)*P_HSF(group,age);
end
    end    

end
for group = 1:4
 errorbar(1:4,P_FF(group,:),prct_err_FF(group,:),'color',colorlist(group),'LineWidth',2);
    hold on;
end
 errorbar(3:4,P_FF(5,3:4),prct_err_FF(5,3:4),'color','r','LineWidth',2);
 
 ylabel(label{1,1});
 set(gca,'Xtick',1:4);
 set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});
   
figure

for group = 1:4
 errorbar(1:4,P_HSF(group,:),prct_err_HSF(group,:),'color',colorlist(group),'LineWidth',2);
    hold on;
end
 errorbar(3:4,P_HSF(5,3:4),prct_err_HSF(5,3:4),'color','r','LineWidth',2);
 
 ylabel(label{1,2});
 set(gca,'Xtick',1:4);
 set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});
    end
%         clear HSF M V M1 V1 x f
        
 








    
