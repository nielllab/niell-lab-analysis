function [meandata se N1 ]=layerAgePlot_frac_responsive(data,ageList,layer,inh,used,used1,label);
ageList=ageList';
colorlist='bmkgrc';

for age=1:4
    for group = 1:5
        if group ==1
            uselist = (ageList==age & (layer<=3) & ~inh & used);
            uselist1=(ageList==age & (layer<=3) & ~inh & used1 );
        elseif group ==2
            uselist = (ageList==age & (layer==4) & ~inh & used);
            uselist1=(ageList==age & (layer==4) & ~inh & used1);
         elseif group==3
             uselist = (ageList==age & (layer==5) & ~inh & used);
             uselist1=(ageList==age & (layer==5) & ~inh & used1);
         elseif group==4
             uselist = (ageList==age & (layer==6) & ~inh &  used);
             uselist1=(ageList==age & (layer==6) & ~inh& used1 );
        elseif group==5
            uselist = (ageList==age & inh & used);
            uselist1=(ageList==age & inh & used1 );
        
        end
        
    if sum(uselist)>=2
        total(group,age)=sum(~isnan(data(uselist1)));
        resp(group,age) = sum(~isnan(data(uselist)));   
        %countdata_lin(group,age) = lin(1,1);
        %N(group,age) = sum(~isnan(data(uselist)));
        frac(group,age)= resp(group,age)./total(group,age);
%       P_FF(group,age)= countdata_FF(group,age)/N(group,age) 
        
        [fr,pci]= binofit(resp(group,age),total(group,age));
         errdata(group,age) = frac(group,age)-pci(1,1);
         sem(group,age)=errdata(group,age)/sqrt(resp(group,age))
%         errdata(group,age) =V/sqrt(total(group,age));
%         prct_err(group,age)= errdata(group,age)/resp(group,age);
%           
%         prct_err_lin(group,age)=prct_err(group,age)*frac(group,age);
       
%            figure
%            hist(data(uselist),0:0.25:2);
    end
    end
end




%Chi_square(2,groups(1,1:2),groups(1,3:4));

frac_E_A=[frac(1:5,1), frac(1:5,4)]
 s_E_A=[sem(1:5,1),sem(1:5,4)]
% 
 figure
 barweb(frac_E_A,s_E_A);

title 'peak1.5hz'

figure
barweb(frac,sem)
ylabel(label);
set(gca,'Xtick',1:5);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','all'});
% legend('EO1','adult');
end










    
