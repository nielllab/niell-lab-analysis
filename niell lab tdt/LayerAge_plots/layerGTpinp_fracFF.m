function [meandata se N1 ]=layerGTpinp_fracFF(data,ageList,layer,inh,pinped,used,label);
ageList=ageList';
colorlist='bmkg';

for age=1:3
  
            uselist = (ageList==age & layer<=5 &  pinped & used);       
        
       
        FF=histc(data(uselist),0:0.005:0.35);
        
        countdata_FF(age)=FF(1,1);
        
       
        N(age) = sum(~isnan(data(uselist)));
        P_FF(age)= countdata_FF(age)/N(age) ;
        
         [fr,pci]= binofit(countdata_FF(age),N(age));
         errdata(age) = P_FF(age)-pci(1,1);
         sem(age)=errdata(age)/sqrt(countdata_FF(age))  ;      

    
end

% figure
% errorbar(1:3,fliplr(frac(1,:)),fliplr(sem(1,:)),'k');hold on
% %errorbar(1:3,frac(1,:),sem(1,:),'k');hold on


%  shadedErrorBar(1:4,frac(group,:),prct_err_lin(group,:),'k');



%Chi_square(2,groups(1,1:2),groups(1,3:4));


figure
barweb(fliplr(P_FF),fliplr(sem))
ylabel(label);
set(gca,'Xtick',1:3);
%set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
legend('WT','N2A','N2B');

% figure
% pie([frac(1) 1-frac(1)])
% title 'N2B'
% figure
% pie([frac(2) 1-frac(2)])
% title 'N2A'
% figure
% pie([frac(3) 1-frac(3)])
% title 'WT'

end
