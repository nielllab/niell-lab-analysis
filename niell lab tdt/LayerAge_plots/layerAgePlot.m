function [meandata errdata N mediandata]=layerAgePlot(data,ageList,layer,inh, used,label);
ageList=ageList';
figure
hold on

colorlist='bgrm';

for age=1:3
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh  & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh  & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh &  used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==5
            uselist = (ageList==age & inh & used);
        elseif group==6
             uselist = (ageList==age & (layer<=6) & ~inh & used);
      
        end
%         clear HSF M V M1 V1 x f
    
        
        
    
        if sum(uselist)>2 
             N(group,age) = sum(~isnan(data(uselist)));
             %meandata(group,age) = nanmean(data(uselist));
             mediandata(group,age)=nanmedian(data(uselist));
             
           
             s = semedian(data(uselist));
             errdata_med(group,age)=s;
             errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
            
        end      
    
    end

end
axis ij
set(gca,'ytick',1:6);
set(gca,'yticklabel',{'2/3','4','5','6','inh','All'});
xlabel(label)
% figure
% barweb(meandata,errdata)
% ylabel(label);
% set(gca,'Xtick',1:5);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
% % 

% Med_E_A=[mediandata(:,1), mediandata(:,4)]
% s_med_E_A=[errdata_med(:,1),errdata_med(:,4)]
% 
% figure
% barweb(Med_E_A,s_med_E_A);



 figure
 barweb(mediandata,errdata_med);
 ylabel(label)
 xlabel('median')
 
%  figure
%  errorbar(1:4,mediandata(1,:),errdata_med(1,:),'k');hold on
%  errorbar(1:4,mediandata(3,:),errdata_med(6,:),'g');hold on
 end
% % % 








    
