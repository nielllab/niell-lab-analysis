function layerAgeScatter(data,ageList,layer,inh,used,label);
ageList=ageList';
figure
hold on
colorlist='bg';
for age=1:2
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh &  used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh & used);
        elseif group==5
            uselist = (ageList==age & inh & used);
        elseif group==6
            uselist = (ageList==age & (layer<=6) & ~inh & used);
        end
        if sum(uselist)>0
            plot(data(uselist),(rand(size(data(uselist)))*0.15 - 0.075)+ones(size(data(uselist)))*(group + (age-1)*0.4 - 0.2),'o','Color',colorlist(age)); hold on
            m  = nanmean(data(uselist));
%             N =  sum(~isnan(data(uselist)));
%             err = nanstd(data(uselist))/sqrt(N);
            
            plot([m m], group + (age-1)*0.4 -0.2 + [ -0.2 0.2],'Color',colorlist(age),'LineWidth',4);
%             plot ([(m-err) (m-err)],group + (age-1)*0.4 -0.2 + [-0.2 0.2],'Color',colorlist(age),'LineWidth',2);
%             plot ([(m+err) (m+err)],group + (age-1)*0.4 -0.2 + [-0.2 0.2],'Color',colorlist(age),'LineWidth',2);
        end
            
    end  
end

axis ij
set(gca,'ytick',1:6);
set(gca,'yticklabel',{'2/3','4','5','6','inh','All'});
xlabel(label)


    
