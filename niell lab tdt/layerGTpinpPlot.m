function [meandata errdata N mediandata]=layerGTpinpPlot(data,pinp,geno,layer,inh,exc,used,label);
pinp=pinp;
geno=geno';


for GT = 1:3
    
for P=1:2
    
    for group = 1:4
        if group ==1
            uselist = (geno==GT & pinp==P & layer<=3 & exc & used);
        elseif group ==2
            uselist = (geno==GT & pinp==P & layer==4 & exc  & used);
        elseif group==3
            uselist = (geno==GT & pinp==P & layer==5 & exc &  used);
        elseif group==4
            uselist = (geno==GT & pinp==P & inh & used);
        end

        if sum(uselist)>1 
             N(group,P) = sum(~isnan(data(uselist)));
             %meandata(group,age) = nanmean(data(uselist));
             mediandata(group,P)=nanmedian(data(uselist));
             s = semedian(data(uselist));
             errdata_med(group,P)=s;
             errdata(group,P)=nanstd(data(uselist))/sqrt(N(group,P));             
             alldata{group,P}= data(uselist);
        else
            mediandata(group,P)=NaN
            errdata_med(group,P)=NaN;
        end
      
    end
end

figure
barweb(mediandata,errdata_med)
ylabel(label);
set(gca,'Xtick',1:4);
set(gca,'Xticklabel',{'2/3','4','5','inh'});
legend('Not-PINPed','PINPed');
title(GT)

end

%  figure
%  errorbar(1:4,mediandata(1,:),errdata_med(1,:),'k');hold on
%  errorbar(1:4,mediandata(3,:),errdata_med(6,:),'g');hold on
 end
% % % 








    
