function [meandata errdata N mediandata]=layerAgePlot(data,ageList,layer,inh, used,label);
ageList=ageList';
figure
hold on

colorlist='bgrm';

for age=1:4
    for group = 1:6
        if group ==1
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh  & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh  & used);
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh &   used);
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &   used);
        elseif group==5
            uselist = (ageList==age & inh  & used);
      
        end
%         clear HSF M V M1 V1 x f
        
        
        %figure
       % hist(data(uselist),0:0.02:.35); %%for use with SF_pref
        %hist(data(uselist),0:0.50:4.5); %%for use with SF_bandwidth
      %  hist(data(uselist),0:1:10)
        
        
   %      [f,x]=hist(data(uselist),0:0.02:.35);%# create histogram from the data
%         %[f,x]= hist(data(uselist),0:0.5:4.5); %%for use with SF_bandwidth
%         %[f,x]=hist(data(uselist),0:45:360);
%           [f,x]=hist(data(uselist),0:0.2:2); %use with F1F0 data
%         % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
%         %figure
%         %bar(x,f/sum(f));
% 
%            
%           countdata_FF(group,age) = histc(data(uselist),0.005); %%to use with SF_pref
%           N(group,age) = sum(~isnan(data(uselist)));
% % % %         
%           P_FF(group,age)= countdata_FF(group,age)/N(group,age)
%           [M,V]= binostat(N(group,age),P_FF(group,age));
%           
%           errdata_FF(group,age) =sqrt(V)/sqrt(N(group,age));
%           prct_err(group,age)= errdata_FF(group,age)/countdata_FF(group,age);
%           prct_err_FF(group,age)=prct_err(group,age)*P_FF(group,age);
% % % %         
%           HSF = histc(data(uselist),0.2:0.15:0.35);
%           countdata_HSF(group,age)= HSF(1,1);
% % % %         
%            P_HSF(group,age)= countdata_HSF(group,age)/N(group,age)
%            [M1,V1]= binostat(N(group,age),P_HSF(group,age));
%            
%            errdata_HSF(group,age) =sqrt(V1)/sqrt(N(group,age));
%            prct_err_H(group,age)= errdata_HSF(group,age)/countdata_HSF(group,age);
%            prct_err_HSF(group,age)=prct_err_H(group,age)*P_HSF(group,age);
% 

        
       
        
    
        if sum(uselist)>2 
             N(group,age) = sum(~isnan(data(uselist)));
             meandata(group,age) = nanmean(data(uselist));
             mediandata(group,age)=nanmedian(data(uselist));
             
           
             s = semedian(data(uselist));
             errdata_med(group,age)=s;
             errdata(group,age)=nanstd(data(uselist))/sqrt(N(group,age));             
             alldata{group,age}= data(uselist);
             
       
            plot(data(uselist),(rand(size(data(uselist)))*0.15 - 0.075)+ones(size(data(uselist)))*(group + (age-1)*0.4 - 0.2),'o','Color',colorlist(age)); hold on
            m  = nanmean(data(uselist));
            s = semedian(data(uselist));
            plot([m m], group + (age-1)*0.4 -0.2 + [ -0.2 0.2],'Color','k','LineWidth',4);
            %plot ([(m-s) (m-s)],group + (age-1)*0.4 -0.2 + [-0.2 0.2],'Color','k','LineWidth',2);
            %plot ([(m+s) (m+s)],group + (age-1)*0.4 -0.2 + [-0.2 0.2],'Color','k','LineWidth',2);
        end
             
    
    end


axis ij
set(gca,'ytick',1:5);
set(gca,'yticklabel',{'2/3','4','5','6','inh','All'});
xlabel(label)  


 
%   figure  
%   barweb(P_FF,prct_err_FF)
% % % % 
%   figure
%   barweb(P_HSF,prct_err_HSF)

end

figure
barweb(meandata,errdata)
ylabel(label);
set(gca,'Xtick',1:5);
set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
legend('EO1','adult');
% % 
 figure
 barweb(mediandata,errdata_med);
 ylabel(label)
 xlabel('median')
 
 end
% % % 


%   ranksum(alldata{1,1},alldata{2,1}) % signfiicance layer2/3
%   ranksum(alldata{3,1},alldata{3,2}) % signfiicance layer4
%   [p h]=kstest2(alldata{6,1},alldata{6,2});
%  ranksum(alldata{3,1},alldata{3,2}) % sign. layer5
%  ranksum(alldata{4,1},alldata{4,2}) % sign layer6
%  ranksum(alldata{5,1},alldata{5,2}) % sign inhibtory
 %ranksum(bothdata{6,1},bothdata{6,2}) % all excitaory

 
 
% else
%     %%do something with bandwidth data??
% end

%end







    
