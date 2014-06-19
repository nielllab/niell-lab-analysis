function [meandata errdata N mediandata]= layer_age_line(data,ageList,layer,inh,used,label );
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ageList=ageList';

colorlist='bmkg';

figure

for group=1:5
    for age =1:4
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
        elseif group==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==5
             uselist = (ageList==age & inh & used);
        end


if sum(uselist)>2 
   

    %N(group,age) = sum(~isnan(data(uselist)));
   N=100
    meandata(group,age) = nanmean(data(uselist));
    errdata(group,age)=nanstd(data(uselist))/sqrt(N);             

    mediandata(group,age)=nanmedian(data(uselist));
    s = semedian(data(uselist));
    errdata_med(group,age)=s;
    alldata{group,age}= data(uselist);
    
end
    end    

end
for group = 1:4
 errorbar(1:4,meandata(group,:),errdata(group,:),'color',colorlist(group),'LineWidth',2);
    hold on;
end
 errorbar(3:4,meandata(5,3:4),errdata(5,3:4),'color','r','LineWidth',2);
 
 ylabel(label);
 set(gca,'Xtick',1:4);
 set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});
   
    
    end
    




