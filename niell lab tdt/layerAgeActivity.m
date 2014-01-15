function [meandata errdata N]=layerAgeActivity(data1,data2,ageList,layer,inh,used,label,titlestr,mid);
ageList=ageList';
for age=1:2
    figure
    
    colors = 'rgbcmy'

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
          plot(data1(uselist),data2(uselist),[colors(group) '.']);hold on
    end
    legend('2/3','4','5','6','inh','mid');
    xlabel(label{1}); ylabel(label{2});
        if age==1
        title(sprintf('%s EO',titlestr));
    else
        title(sprintf('%s adult',titlestr));
        end
   axis square
   axis equal
    
    plot([0 10],[0 10])
end



    
