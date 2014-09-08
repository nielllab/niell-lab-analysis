function [P_FF P_HSF prct_err_FF prct_err_HSF data N ]=layerAgePlot_SF_line(data,ageList,layer,inh,used,label);
ageList=ageList';
figure
hold on

colorlist='bmkgr';

figure

for group=1:3
    for age =1:4
        if group ==1
            uselist = (ageList==age & (layer<=4) & ~inh & used);
        elseif group==2
            uselist = (ageList==age & (layer>=5) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer<=6) & ~inh & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh &  used);
%         elseif group==3
%              uselist = (ageList==age & inh & used);
        end


if sum(uselist)>2 
   

    [f,x]=hist(data(uselist),0:0.02:.35);%# create histogram from the data
          
          countdata_FF(group,age) = histc(data(uselist),0.005); %%to use with SF_pref
          N(group,age) = sum(~isnan(data(uselist)));
% % % %         
          P_FF(group,age)= countdata_FF(group,age)/N(group,age)
         
          %[M,V]= binostat(N(group,age),P_FF(group,age));
          [fr,pci]= binofit(countdata_FF(group,age),N(group,age));
         errdata(group,age) = P_FF(group,age)-pci(1,1);
         sem_FF(group,age)=errdata(group,age)/sqrt(N(group,age))
          
%           errdata_FF(group,age) =sqrt(V)/sqrt(N(group,age));
%           prct_err(group,age)= errdata_FF(group,age)/countdata_FF(group,age);
%           
%           prct_err_FF(group,age)=prct_err(group,age)*P_FF(group,age);
% % % %         
          HSF = histc(data(uselist),0.3:0.05:0.35);
          countdata_HSF(group,age)= HSF(1,1);
% % % %         
          P_HSF(group,age)= countdata_HSF(group,age)/N(group,age)
       
          [fr,pci_H]= binofit(countdata_HSF(group,age),N(group,age));
         errdata_H(group,age) = P_HSF(group,age)-pci_H(1,1);
         sem_H(group,age)=errdata_H(group,age)/sqrt(N(group,age))
end
    end    

end
% for group = 1:3
%  errorbar(1:4,P_FF(group,:),prct_err_FF(group,:),'color',colorlist(group),'LineWidth',2);
%     hold on;
% end

figure 
errorbar(1:4,P_FF(1,:),sem_FF(1,:),'k');hold on
errorbar(1:4,P_FF(2,:),sem_FF(2,:),'g');hold on
% shadedErrorBar(1:4,P_FF(1,:),prct_err_FF(1,:),'k');hold on
% shadedErrorBar(1:4,P_FF(2,:),prct_err_FF(2,:),'b');hold on
% 
% shadedErrorBar(3:4,P_FF(3,3:4),prct_err_FF(3,3:4),'r');
 

   
figure
 errorbar(1:4,P_HSF(3,:),sem_H(3,:),'k');hold on
 %errorbar(1:4,P_HSF(2,:),sem_H(2,:),'g');hold on
% for group = 1:4
%  errorbar(1:4,P_HSF(group,:),prct_err_HSF(group,:),'color',colorlist(group),'LineWidth',2);
%     hold on;
% end
% shadedErrorBar(1:4,P_HSF(1,:),prct_err_HSF(1,:),'k');hold on 
% shadedErrorBar(1:4,P_HSF(2,:),prct_err_HSF(2,:),'b');hold on
% shadedErrorBar(3:4,P_HSF(3,3:4),prct_err_HSF(3,3:4),'r'); 
 
    end
%         clear HSF M V M1 V1 x f
        
 








    
