function [meandata se N1 ]=layerAgePlot_frac_responsive(data,ageList,layer,inh,used,label);
ageList=ageList';
colorlist='bmkgrc';

for age=1:2
    for group = 1:5
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
            uselist1=(ageList==age & (layer<=3) & ~inh );
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
            uselist1=(ageList==age & (layer==4) & ~inh );
        elseif group==3
            uselist = (ageList==age & (layer==5) & ~inh & used);
            uselist1=(ageList==age & (layer==5) & ~inh );
        elseif group==4
            uselist = (ageList==age & (layer==6) & ~inh &  used);
            uselist1=(ageList==age & (layer==6) & ~inh );
        elseif group==5
            uselist = (ageList==age & inh & used);
            uselist1=(ageList==age & inh );
      
        end
        
    if sum(uselist)>=2
        total(group,age)=sum(~isnan(data(uselist1)));
        resp(group,age) = sum(~isnan(data(uselist)));   
        %countdata_lin(group,age) = lin(1,1);
        %N(group,age) = sum(~isnan(data(uselist)));
        frac(group,age)= resp(group,age)/total(group,age);
%       P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
        [M,V]= binostat(total(group,age),frac(group,age));
          
        errdata(group,age) =V/sqrt(total(group,age));
        prct_err(group,age)= errdata(group,age)/resp(group,age);
          
        prct_err_lin(group,age)=prct_err(group,age)*frac(group,age);
       
%            figure
%            hist(data(uselist),0:0.25:2);
    end
    end
end

figure

for group = 1:4
 errorbar(1:2,frac(group,:),prct_err_lin(group,:),'color',colorlist(group),'LineWidth',2);
    hold on;
end
 
%errorbar(3:4,frac(5,3:4),prct_err_lin(5,3:4),'color','r','LineWidth',2);

%  ylabel(label{1,1});
%  set(gca,'Xtick',1:4);
%  set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});

% ranksum(bothdata{1,1},bothdata{1,2}) % signfiicance layer2/3
% ranksum(bothdata{2,1},bothdata{2,2}) % signfiicance layer4
% ranksum(bothdata{3,1},bothdata{3,2}) % sign. layer5
% ranksum(bothdata{4,1},bothdata{4,2}) % sign layer6
% % ranksum(bothdata{5,1},bothdata{5,2})% sign inhibotyr
% ranksum(bothdata{6,1},bothdata{6,2}) % all excitaor

% groups =zeros(6,4);
% groups(:,1)= bothdata_lin(:,1);
% groups(:,2)=  bothdata_other(:,1);
% groups(:,3)=  bothdata_lin(:,2);
% groups(:,4) = bothdata_other(:,2);
% Group1 = groups(1,1:2)
% Group2 = groups(1,3:4)


%Chi_square(2,groups(1,1:2),groups(1,3:4));


figure
barweb(frac,prct_err_lin)
ylabel(label);
set(gca,'Xtick',1:5);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
end










    
