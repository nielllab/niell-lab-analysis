function [mediandata errdata ratio]  =layerAgePlot_ratio_jlh(data1,data2,ageList,layer,inh,used,label,titlestr)
%plot ratio of data1 vs data2 with CIs by layer

ageList=ageList';
 for age=1:2
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
        elseif group==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh & used);
        elseif group==5   
             uselist = (ageList==age & inh & used);
        elseif group==6
            uselist = (ageList==age & (layer<=6) & ~inh & used);
        end
   
        if sum(uselist)>2
        
        ratio_D2_D1(group,age)=nanmedian(data2(uselist))/nanmedian(data1(uselist));
       
        sem_ratio=semedian_ratio(data2(uselist),data1(uselist));
        err_ratio(group,age)=sem_ratio;
        
        %rat(group, age)=nanmedian(data2(uselist)./data1(uselist));
        mx=nanmedian(data1(uselist));
        sx= semedian(data1(uselist));
                     
        my=nanmedian(data2(uselist));
        sy=semedian(data2(uselist));
%                      
        mediandata_x(group,age)=mx;
        errdata_med_x(group,age)=sx;
                     
        mediandata_y(group,age)=my;
        errdata_med_y(group,age)=sy;
        
        end
    end
 end
 
 
 
 
Med_E_A=[ratio_D2_D1(:,1), ratio_D2_D1(:,2)]
s_med_E_A=[err_ratio(:,1),err_ratio(:,2)]

figure
barweb(Med_E_A,s_med_E_A);
 
figure
barweb(ratio_D2_D1,err_ratio)



% figure
% barweb(rat,err_ratio)

 figure
 errorbar(1:4,ratio_D2_D1(1,:),err_ratio(1,:),'k');hold on
 errorbar(1:4,ratio_D2_D1(2,:),err_ratio(2,:),'r');hold on
%  
figure
errorbar(1:4,mediandata_x(5,:),errdata_med_x(5,:),'k');hold on 
errorbar(1:4,mediandata_y(5,:), errdata_med_y(5,:),'g');hold on
title ' layer 5 spont running_green stat_black'

% figure
% errorbar(1:4,mediandata_x(1,:),errdata_med_x(1,:),'k');hold on 
% errorbar(1:4,mediandata_y(1,:), errdata_med_y(1,:),'g');hold on
% title ' superficial'
end

