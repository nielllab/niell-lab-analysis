function [meandata se N1 ]=layerAgePlot_frac_simple(data,ageList,layer,inh,used,label);
ageList=ageList';
for age=1:2
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==5
            uselist = (ageList==age & inh & used);
        elseif group==6
            uselist = (ageList==age & (layer<=6)& ~inh & used);
        end
        clear h V phat pci lin N frac M other
        
        h = data(uselist);
%             x = 0.0:0.2:2.0;
%             hist(h,x);
%             axis xy;
%             xlabel(label);
%             ylabel(age);
%             title 'layer 2/3';
           N = sum(~isnan(data(uselist))); 
           lin = histc(h,1:1:2);
           %other = N-lin;
           frac = (lin(1,1)/N);
           
          % frac_lin(group,age) = (lin(1,1)/totalct(1,1));
           
         [phat,pci]= binofit(lin(1,1),N);
         meandata(group,age) = phat
         N1(group,age) = sum(~isnan(data(uselist)));
         se(group,age) = sqrt((phat*(1-phat))/N); %%working out standard error for binomial and multinomial distributions
         %bothdata_lin(group,age)= frac
         %bothdata_other(group,age)= 1-frac;
         %CI(group,age) = pci;
         %bothdata{group,age}= data(uselist);
            
    end
end
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
barweb(meandata,se)
ylabel(label);
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
legend('EO1','adult');










    
