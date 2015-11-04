function [mediandata errdata ratio]  =layerGenoPlot_pinp_ratio(data1,data2,pinp,geno,layer,used,label)
%plot ratio of data1 vs data2 with CIs by layer
pinp=pinp;
geno=geno';

for GT = 1:3
    
for P=1:2
    
    for group = 1:2
        if group ==1
            uselist = (geno==GT & pinp==P & layer<=4 & used);
        elseif group==2
            uselist = (geno==GT & pinp==P  & used);
%         elseif group==3
%             uselist = (geno==GT & pinp==P & layer==5 & exc & used);
%         elseif group==4  
%              uselist = (geno==GT & pinp==P & inh & used);
        end
   
        if sum(uselist)>1
        
        ratio_D2_D1(group,P)=nanmedian(data2(uselist))/nanmedian(data1(uselist));
        sem_ratio=semedian_ratio(data2(uselist),data1(uselist));
        err_ratio(group,P)=sem_ratio;
        N(group,P)=sum(~isnan(data1(uselist)));
        mx=nanmedian(data1(uselist));
        sx= semedian(data1(uselist));
                     
        my=nanmedian(data2(uselist));
        sy=semedian(data2(uselist));
%                      
        mediandata_x(group,P)=mx;
        errdata_med_x(group,P)=sx;
                     
        mediandata_y(group,P)=my;
        errdata_med_y(group,P)=sy;
        
        else 
          ratio_D2_D1(group,P)=NaN;
          err_ratio(group,P)=NaN;
          mediandata_x(group,P)=NaN;
          errdata_med_x(group,P)=NaN;
          mediandata_y(group,P)=NaN;
          errdata_med_y(group,P)=NaN;
        end
    end
    
 end
 
figure
barweb(ratio_D2_D1,err_ratio)
ylabel(label);
set(gca,'Xtick',1:4);
set(gca,'Xticklabel',{'2/3','4','5','inh'});
legend('Not-PINPed','PINPed');
title(GT)

end 
% Med_E_A=[ratio_D2_D1(:,1), ratio_D2_D1(:,2)]
% s_med_E_A=[err_ratio(:,1),err_ratio(:,2)]
% 
% figure
% barweb(Med_E_A,s_med_E_A);
 
% 
%  figure
%  errorbar(1:4,ratio_D2_D1(1,:),err_ratio(1,:),'k');hold on
%  errorbar(1:4,ratio_D2_D1(2,:),err_ratio(2,:),'r');hold on
% %  
% figure
% errorbar(1:4,mediandata_x(5,:),errdata_med_x(5,:),'k');hold on 
% errorbar(1:4,mediandata_y(5,:), errdata_med_y(5,:),'g');hold on
% title ' layer 5 spont running_green stat_black'

% figure
% errorbar(1:4,mediandata_x(1,:),errdata_med_x(1,:),'k');hold on 
% errorbar(1:4,mediandata_y(1,:), errdata_med_y(1,:),'g');hold on
% title ' superficial'
end

