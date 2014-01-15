function [meandata errdata N]=layerAgePlotMv(data,ageList,layer,inh,used,label,mid);
ageList=ageList';
for mv = 1:2
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
            meandata(group,(age-1)*2 + mv) = nanmedianMW(data(uselist,mv));
            N(group,(age-1)*2 + mv) = sum(~isnan(data(uselist,mv)));
        errdata(group,(age-1)*2 + mv) = nanstd(data(uselist,mv))/sqrt(N(group,(age-1)*2+mv));
    end
    end
end

figure
barweb(meandata,errdata);
ylabel(label);
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',{'2/3','4','5','6','inh','mid'});
legend('EO1 -stop','EO1 - mv','adult -stop','adult - mv');


    
