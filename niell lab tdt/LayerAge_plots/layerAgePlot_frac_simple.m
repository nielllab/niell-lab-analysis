function [meandata se N1 ]=layerAgePlot_frac_simple(data,ageList,layer,inh,pinped,used,label);
ageList=ageList';
colorlist='bmkg';

for age=1:3
    %for group = 1:3
      %  if group ==1
            uselist = (ageList==age & layer<=5 &  pinped & used);
%         elseif group ==2
%             uselist = (ageList==age & (layer==4|layer==3) & ~inh & pinped & used);
%         elseif group==3
%             uselist = (ageList==age & (layer==5) & ~inh & pinped & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh & pinped &  used);
%         elseif group==5
%             uselist = (ageList==age & inh & used);
%         elseif group==6
%             uselist = (ageList==age & (layer<=6)& ~inh & used);
      %  end
        
        lin = histc(data(uselist),0.8:1:2); 
        countdata_lin(age) = lin(1,1);
       
        N(age) = sum(~isnan(data(uselist)));
        frac(age)= countdata_lin(age)/N(age);
        
        [fr,pci]= binofit(countdata_lin(age),N(age));

         errdata(age) = frac(age)-pci(1,1);

        sem(age)=errdata(age)/sqrt(countdata_lin(age))        

       
%            figure
%            hist(data(uselist),0:0.25:2);
            
    %end
end

% figure
% errorbar(1:3,fliplr(frac(1,:)),fliplr(sem(1,:)),'k');hold on
% %errorbar(1:3,frac(1,:),sem(1,:),'k');hold on



%Chi_square(2,groups(1,1:2),groups(1,3:4));


figure
barweb(fliplr(frac),fliplr(sem))
ylabel(label);
set(gca,'Xtick',1:3);
%set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
legend('WT','N2A','N2B');



figure
pie([frac(1) 1-frac(1)])
title 'N2B'
figure
pie([frac(2) 1-frac(2)])
title 'N2A'
figure
pie([frac(3) 1-frac(3)])
title 'WT'


end










    
