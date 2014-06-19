function [meandata errdata N mediandata]=layerAgePlot(data,ageList,layer,inh,used,label,titlestr);
ageList=ageList';
for age=1:4
    for group = 1:2
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & used);
%         elseif group ==2
%             uselist = (ageList==age & (layer==4) & ~inh & used);
%         elseif group==3
%             uselist = (ageList==age & (layer==5) & ~inh & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh &  used);
         elseif group==2
             uselist = (ageList==age & inh & used);
%        elseif group==6
%            uselist = (ageList==age & (layer<=6)& ~inh & used);
        end
        clear HSF M V M1 V1 x f
        
        
        figure
        hist(data(uselist),0:45:330)
        N(group,age) = sum(~isnan(data(uselist)));
        xlabel(label{1}); ylabel(label{2});
        if age==1
        title(sprintf('%s EO',titlestr));
        elseif age==2
        title(sprintf('%s EO3',titlestr));
        elseif age==3
        title(sprintf('%s EO7',titlestr));
        elseif age==4
        title(sprintf('%s adult',titlestr));
        end
        
        [f,x]=hist(data(uselist),0:45:330);
        % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
        
        figure
        bar(x,f/sum(f));
         xlabel(label{1}); ylabel(label{2});
        if age==1
        title(sprintf('%s EO',titlestr));
        elseif age==2
        title(sprintf('%s EO3',titlestr));
        elseif age==3
        title(sprintf('%s EO7',titlestr));
        elseif age==4
        title(sprintf('%s adult',titlestr));
        end
     

%               meandata(group,age) = nanmean(data(uselist));
%               errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
%               mediandata(group,age) = nanmedian(data(uselist));
% %              errdataMed(group,age) = mad(data(uselist),1);
%              bothdata{group,age}= data(uselist);
             
        end
    end

   



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







    
