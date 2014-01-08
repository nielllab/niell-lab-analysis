function [meandata errdata N mediandata]=layerAgePlot(data,ageList,layer,inh,used,label);
ageList=ageList';
for age=1:2
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
        elseif group==5
            uselist = (ageList==age & inh & used);
        elseif group==6
            uselist = (ageList==age & (layer<=6)& ~inh & used);
        end
        clear HSF M V M1 V1 x f
        
        
        figure
        %hist(data(uselist),0:0.02:.32); %%for use with SF_pref
        %hist(data(uselist),0:0.50:4.5); %%for use with SF_bandwidth
        hist(data(uselist),0:45:360)
        N(group,age) = sum(~isnan(data(uselist)));
        
        
        %[f,x]=hist(data(uselist),0:0.02:.35);%# create histogram from the data
        %[f,x]= hist(data(uselist),0:0.5:4.5); %%for use with SF_bandwidth
        [f,x]=hist(data(uselist),0:45:360);
        % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
        figure
        bar(x,f/sum(f));
        
        %if exist ('wpref_stat', 'var')
           
        %countdata_FF(group,age) = histc(data(uselist),0.005); %%to use with SF_pref
        
%         countdata_FF(group,age) = histc(data(uselist),4);
%         
%         P_FF(group,age)= countdata_FF(group,age)/N(group,age)
%         [M,V]= binostat(N(group,age),P_FF(group,age));
%         errdata_FF(group,age) =sqrt(V)/sqrt(N(group,age));
%         prct_err(group,age)= errdata_FF(group,age)/countdata_FF(group,age);
%         prct_err_FF(group,age)=prct_err(group,age)*P_FF(group,age);
%         
%         %HSF = histc(data(uselist),0.24:0.11:0.35);
%         HSF = histc(data(uselist),4.5);
%         countdata_HSF(group,age)= HSF(1,1);
%         
%         P_HSF(group,age)= countdata_HSF(group,age)/N(group,age)
%         [M1,V1]= binostat(N(group,age),P_HSF(group,age));
%         errdata_HSF(group,age) =sqrt(V1)/sqrt(N(group,age));
%         prct_err_H(group,age)= errdata_HSF(group,age)/countdata_HSF(group,age);
%         prct_err_HSF(group,age)=prct_err_H(group,age)*P_HSF(group,age);

        %elseif exist ('driftwbw_1', 'var')
        
       
        
        %else
%              meandata(group,age) = nanmean(data(uselist));
%              errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
%              mediandata(group,age) = nanmedian(data(uselist));
%              errdataMed(group,age) = mad(data(uselist),1);
%              bothdata{group,age}= data(uselist);
             
        end
    end

   


%if exist ('P_FF','var')
    
% figure  
% barweb(P_FF,prct_err_FF)
% 
% figure
% barweb(P_HSF,prct_err_HSF)

% elseif exist('meandata','var')
% 
% figure
% barweb(meandata,errdata)
% ylabel(label);
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
% % 
% figure
% barweb(mediandata,errdataMed);
% ylabel(label)
% xlabel('median')
% % % 
%  ranksum(bothdata{1,1},bothdata{1,2}) % signfiicance layer2/3
%  ranksum(bothdata{2,1},bothdata{2,2}) % signfiicance layer4
%  ranksum(bothdata{3,1},bothdata{3,2}) % sign. layer5
%  ranksum(bothdata{4,1},bothdata{4,2}) % sign layer6
%  ranksum(bothdata{5,1},bothdata{5,2}) % sign inhibtory
%  ranksum(bothdata{6,1},bothdata{6,2}) % all excitaory

% else
%     %%do something with bandwidth data??
% end

end







    
