function layerAgeScatter(data,ageList,layer,inh,used,label,mid);
ageList=ageList';
figure
hold on
colorlist='br';
for age=1:2
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & ~mid & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & ~mid& used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & ~mid& used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh & ~mid& used);
        elseif group==5
            uselist = (ageList==age &inh & used);
        elseif group==6
            uselist = (ageList==age &mid & used);
        end
        if sum(uselist)>0
            plot(data(uselist),ones(size(data(uselist)))*(group + (age-1)*0.3),'o','Color',colorlist(age))
        end
    end
end

axis ij
set(gca,'ytick',1:6);
set(gca,'yticklabel',{'2/3','4','5','6','inh','mid'});
xlabel(label)


    
