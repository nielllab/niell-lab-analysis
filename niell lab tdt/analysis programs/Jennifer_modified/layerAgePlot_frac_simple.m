function [meandata se N1 ]=layerAgePlot_frac_simple(data,ageList,layer,inh,pinped, used,label);
ageList=ageList';
colorlist='bmkg';

for age=1:3
    for group = 1:2
        if group ==1
            uselist = (ageList==age & (layer<=5) & ~inh & pinped & used);
        elseif group ==2
            uselist = (ageList==age & (layer<=5) & ~inh &~pinped & used);
%         elseif group==3
%             uselist = (ageList==age & (layer==5) & ~inh & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh &  used);
%         elseif group==5
%             uselist = (ageList==age & inh & used);
%         elseif group==6
%             uselist = (ageList==age & (layer<=4)& ~inh & used);
        end
        
        
        
        lin = histc(data(uselist),0.8:1:2);   
        countdata_lin(group,age) = lin(1,1);
        N(group,age) = sum(~isnan(data(uselist)));
        frac(group,age)= countdata_lin(group,age)/N(group,age);
%           P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
        [fr,pci]= binofit(countdata_lin(group,age),N(group,age));
         errdata(group,age) = frac(group,age)-pci(1,1);
         sem(group,age)=errdata(group,age)/sqrt(countdata_lin(group,age))        

       
%            figure
%            hist(data(uselist),0:0.25:2);
            
    end
end

% figure
% errorbar(1:3,frac(1,:),sem(1,:),'k');hold on
% errorbar(1:3,frac(2,:),sem(2,:),'g');hold on
% 
% figure
% errorbar(1:3,frac(6,:),sem(6,:),'k');hold on
% errorbar(1:3,frac(1,:),sem(1,:),'k');hold on
% for group = 1:2
%  shadedErrorBar(1:4,frac(group,:),prct_err_lin(group,:),'k');
%     hold on;
% end
 
%  ylabel(label{1,1});
%  set(gca,'Xtick',1:4);
%  set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});

% ranksum(bothdata{1,1},bothdata{1,2}) % signfiicance layer2/3
% ranksum(bothdata{2,1},bothdata{2,2}) % signfiicance layer4
% ranksum(bothdata{3,1},bothdata{3,2}) % sign. layer5
% ranksum(bothdata{4,1},bothdata{4,2}) % sign layer6
% % ranksum(bothdata{5,1},bothdata{5,2})% sign inhibotyr
% ranksum(bothdata{6,1},bothdata{6,2}) % all excitaor

% groups =zeros(6,4);
% groups(:,1)= bothdata_lin(:,1);
% groups(:,2)=  bothdata_other(:,1);
% groups(:,3)=  bothdata_lin(:,2);
% groups(:,4) = bothdata_other(:,2);
% Group1 = groups(1,1:2)
% Group2 = groups(1,3:4)


%Chi_square(2,groups(1,1:2),groups(1,3:4));


figure
barweb(frac,sem)
ylabel(label);
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',{'pinped','non_pinped'});
legend('N2B','N2A','wt');
end










    
