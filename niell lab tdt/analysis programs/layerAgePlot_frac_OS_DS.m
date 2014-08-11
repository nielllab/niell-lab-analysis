function [meandata se N1 ]=layerAgePlot_frac_OS_DS(data,ageList,layer,inh,used,used1,label);
ageList=ageList';
colorlist='bmkgrc';

for age=1:4
    for group = 1:5
        if group ==1
            uselist = (ageList==age & (layer<=3 ) & ~inh & used);
            uselist1=(ageList==age & (layer<=3) & ~inh & used1);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
            uselist1=(ageList==age & (layer==4) & ~inh & used1);

        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
            uselist1=(ageList==age & (layer==5) & ~inh & used1);

        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
            uselist1=(ageList==age & (layer==6) & ~inh & used1 );

        elseif group==5
            uselist = (ageList==age & inh & used);
            uselist1=(ageList==age & inh & used1 );

      
        end
        
        
        total(group,age)=sum(~isnan(data(uselist)));
        resp(group,age) = sum(~isnan(data(uselist1)));   
        %countdata_lin(group,age) = lin(1,1);
        %N(group,age) = sum(~isnan(data(uselist)));
        frac(group,age)= resp(group,age)/total(group,age);
%           P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
         [fr,pci]= binofit(resp(group,age),total(group,age));
         errdata(group,age) = frac(group,age)-pci(1,1);
         sem(group,age)=errdata(group,age)/sqrt(resp(group,age))
            
    end
end

figure
errorbar(1:4,frac(1,:),sem(1,:),'b');hold on
errorbar(1:4,frac(2,:),sem(2,:),'k');hold on

% shadedErrorBar(1:4,frac(3,:),sem(3,:),'g');hold on


% for group = 1:4
%  errorbar(1:4,frac(group,:),sem(group,:),'color',colorlist(group),'LineWidth',2);
%     hold on;
% end
 
% errorbar(3:4,frac(5,3:4),sem(5,3:4),'color','r','LineWidth',2);

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
set(gca,'Xtick',1:5);
set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
legend('EO1','adult');
end










    
