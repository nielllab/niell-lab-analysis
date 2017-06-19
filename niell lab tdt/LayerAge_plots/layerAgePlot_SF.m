function [meandata errdata N mediandata]=layerAgePlot_SF(data,ageList,layer,inh,pinp,used,label);
ageList=ageList';
figure
hold on

colorlist='bmrg';

for age=1:3
   
         uselist = (ageList==age & pinp & used);
     
         clear HSF M V M1 V1 x f
        
        
     
           
          countdata_FF(age) = histc(data(uselist),0.005); %%to use with SF_pref
          N(age) = sum(~isnan(data(uselist))); 
          frac_FF(age)= countdata_FF(age)/N(age);
          
          [fr,pci]= binofit(countdata_FF(age),N(age)); 
          errdata(age) =frac_FF(age)-pci(1,1);
          sem(age)=errdata(age)/sqrt(N(age));
          
          
% %         
          HSF = histc(data(uselist),0.18:0.18:0.36);
          countdata_HSF(age)= HSF(1,1);  
          frac_HSF(age)= countdata_HSF(age)/N(age);
           
          [fr_H,pci_H]= binofit(countdata_HSF(age),N(age)); 
          errdata_H(age) =frac_HSF(age)-pci_H(1,1);
          sem_H(age)=errdata_H(age)/sqrt(N(age));

%         [fr,pci]= binofit(LP(group,age),N(group,age)); 
%         errdata(group,age) =frac(group,age)-pci(1,1);
%         sem(group,age)=errdata(group,age)/sqrt(N(group,age))
       

end

figure
barweb(fliplr(frac_FF),fliplr(sem))
ylabel(label);
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
% % 
 figure
 barweb(fliplr(frac_HSF),fliplr(sem_H));
 ylabel(label)
 xlabel('median')
 
 end
% % % 


%   ranksum(alldata{1,1},alldata{2,1}) % signfiicance layer2/3
%   ranksum(alldata{3,1},alldata{3,2}) % signfiicance layer4
%   [p h]=kstest2(alldata{6,1},alldata{6,2});
%  ranksum(alldata{3,1},alldata{3,2}) % sign. layer5
%  ranksum(alldata{4,1},alldata{4,2}) % sign layer6
%  ranksum(alldata{5,1},alldata{5,2}) % sign inhibtory
 %ranksum(bothdata{6,1},bothdata{6,2}) % all excitaory

 
 
% else
%     %%do something with bandwidth data??
% end

%end







    
