function [meandata errdata N mediandata]= layer_age_line(data,ageList,layer,inh,used,label );
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ageList=ageList';

colorlist='bmkg';

figure

for group=1:4
    for age =1:4
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
         elseif group==2
             uselist = (ageList==age & (layer==4) & ~inh & used);
      elseif group==3
             uselist = (ageList==age & (layer==5) & ~inh & used);
% %         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==4
             uselist = (ageList==age & inh & used);
        end


if sum(uselist)>2 
   

    N(group,age) = sum(~isnan(data(uselist)));
%     N=100
%     meandata(group,age) = nanmean(data(uselist));
%     errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             

    mediandata(group,age)=nanmedian(data(uselist));
    s = semedian(data(uselist));
    errdata_med(group,age)=s;
    alldata{group,age}= data(uselist);
    
end
    end    

end
% for group = 1:3
%  shadedErrorBar(1:4,meandata(group,:),errdata(group,:),'k');
%     hold on;
% end

figure
errorbar(1:4,mediandata(1,:),errdata_med(1,:),'K');hold on
%errorbar(1:4,mediandata(2,:),errdata_med(2,:),'g');hold on
errorbar(1:4,mediandata(3,:),errdata_med(3,:),'b');hold on
errorbar(3:4,mediandata(4,3:4),errdata_med(4,3:4),'r');hold on


%  shadedErrorBar(1:4,meandata(2,:),errdata(2,:),'k');hold on

%shadedErrorBar(1:4,mediandata(2,:),errdata_med(2,:),'k');hold on

%shadedErrorBar(3:4,meandata(3,3:4),errdata(3,3:4),'r');
 
 ylabel(label);
 set(gca,'Xtick',1:4);
 set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});
% for group = 1:4
%  errorbar(1:4,meandata(group,:),errdata(group,:),'color',colorlist(group));
%     hold on;
% end
%  errorbar(3:4,meandata(5,3:4),errdata(5,3:4),'r');
  
    
    end
    




