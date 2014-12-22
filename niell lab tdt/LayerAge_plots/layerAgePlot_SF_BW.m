function [meandata errdata N mediandata]=layerAgePlot_SF_BW(data,ageList,layer,inh,used,label,titlestr);
ageList=ageList';
for age=1:4
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer<=4) & ~inh & used);
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
        clear HSF M V M1 V1 x f
        
        
%         if sum(uselist)>2
        N(group,age) = sum(~isnan(data(uselist)));
        [f,x]=hist(data(uselist),0:1:10);
        LP(group,age)=f(1,10);
        frac(group,age)= LP(group,age)/N(group,age);

        
        
        % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
        
      
%         [fr,pci]= binofit(resp(group,age),total(group,age));
%          errdata(group,age) = frac(group,age)-pci(1,1);
%          sem(group,age)=errdata(group,age)/sqrt(resp(group,age))
        
        [fr,pci]= binofit(LP(group,age),N(group,age)); 
        errdata(group,age) =frac(group,age)-pci(1,1);
        sem(group,age)=errdata(group,age)/sqrt(N(group,age))
        
        [f1,x1]=hist(data(uselist),0:1:10);
        HP(group,age)=f1(1,11);
        frac_H(group,age)= HP(group,age)/N(group,age);

        [fr_1,pci_1]= binofit(HP(group,age),N(group,age)); 
        errdata_H(group,age) =frac_H(group,age)-pci_1(1,1);
        sem_H(group,age)=errdata_H(group,age)/sqrt(N(group,age))

%               meandata(group,age) = nanmean(data(uselist));
%               errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
%               mediandata(group,age) = nanmedian(data(uselist));
% %              errdataMed(group,age) = mad(data(uselist),1);
%              bothdata{group,age}= data(uselist);
%         end    
        end
    end

   
figure
barweb(frac,sem)

figure
barweb(frac,errdata)

figure
errorbar(1:4,frac(1,:),sem(1,:),'K');hold on

figure
barweb(frac_H,sem_H)

% barweb(meandata,errdata)
% ylabel(label);
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
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







    
