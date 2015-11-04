function [meandata errdata N mediandata]=layerAgePlot_pref_Orient(data,ageList,layer,inh,pinped,used,label,titlestr);
ageList=ageList';
for age=1:3
    for group = 1:2
        if group ==1
            uselist = (ageList==age  & pinped & used);
         elseif group ==2
           uselist = (ageList==age  & ~pinped & used);
%         elseif group==3
%             uselist = (ageList==age & (layer==5) & ~inh & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh &  used);
%          elseif group==2
%              uselist = (ageList==age & inh & used);
%        elseif group==6
%            uselist = (ageList==age & (layer<=6)& ~inh & used);
        end
        clear HSF M V M1 V1 x f
        
        
        figure
        hist(data(uselist),0:45:330)
        N(group,age) = sum(~isnan(data(uselist)));
%         xlabel(label{1}); ylabel(label{2});
%         if age==1
%         title(sprintf('%s EO',titlestr));
%         elseif age==2
%         title(sprintf('%s EO3',titlestr));
%         elseif age==3
%         title(sprintf('%s EO7',titlestr));
%         elseif age==4
%         title(sprintf('%s adult',titlestr));
%         end
        
        [f,x]=hist(data(uselist),0:45:330);
        % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
        card = [1 3 5 7];
        s_card(group,age)=sum(f(card));
        frac_card(group,age)= s_card(group,age)/N(group,age);
%           P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
        [fr,pci]= binofit(s_card(group,age),N(group,age));
         errdata(group,age) = frac_card(group,age)-pci(1,1);
         sem(group,age)=errdata(group,age)/sqrt(N(group,age))

        
          
       
     

%               meandata(group,age) = nanmean(data(uselist));
%               errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
%               mediandata(group,age) = nanmedian(data(uselist));
% %              errdataMed(group,age) = mad(data(uselist),1);
%              bothdata{group,age}= data(uselist);
             
        end
    end

   
% figure
% errorbar(1:3,frac_card(1,:),sem(1,:),'k');hold on

%errorbar(1:4,frac_card(2,:),prct_err_lin(2,:),'k');hold on

figure
barweb(frac_card(:,:),sem(:,:))
ylabel(label);
set(gca,'Xtick',1:2);
set(gca,'Xticklabel',{'pinped','non_pinped'});
legend('wt','N2A','N2B');
% % 
% figure
% barweb(mediandata,errdataMed);
% ylabel(label)
% xlabel('median')
% % % 
%  ranksum(bothdata{1,1},bothdata{1,2}) % signfiicance layer2/3
%  ranksum(bothdata{2,1},bothdata{2,2}) % signfiicance layer4
%  ranksum(bothdata{3,1},bothdata{3,2}) % sign. layer5
%  ranksum(bothdata{4,1},bothdata{4,2}) % sign layer6
%  ranksum(bothdata{5,1},bothdata{5,2}) % sign inhibtory
%  ranksum(bothdata{6,1},bothdata{6,2}) % all excitaory

% else
%     %%do something with bandwidth data??
% end

end







    
