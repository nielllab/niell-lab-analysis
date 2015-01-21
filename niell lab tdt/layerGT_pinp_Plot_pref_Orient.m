function [meandata errdata N mediandata]=layerGT_pinp_Plot_pref_Orient(data,pinp,geno,layer,inh,exc,used,label);
pinp=pinp;
geno=geno';

for GT=1:3

for P=1:2   
    for group = 1:3
        if group ==1
            uselist = (geno==GT & pinp==P & layer<=3 & exc & used);
        elseif group==2
            uselist = (geno==GT & pinp==P & layer==4 & exc & used);
        elseif group==3
            uselist = (geno==GT & pinp==P & layer==5 & exc & used);
        
        end
   
        if sum(uselist)>1
       
        figure
        hist(data(uselist),0:45:330)
        N(group,P) = sum(~isnan(data(uselist)));
        title (GT)
        [f,x]=hist(data(uselist),0:45:330);
        % g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
        card = [1 3 5 7];
        s_card(group,P)=sum(f(card));
        frac_card(group,P)= s_card(group,P)/N(group,P);
%           P_FF(group,P)= countdata_FF(group,P)/N(group,P) 
        [fr,pci]= binofit(s_card(group,P),N(group,P));
        errdata(group,P) = frac_card(group,P)-pci(1,1);
        sem(group,P)=errdata(group,P)/sqrt(N(group,P));
             
        end
    end

end
% figure
% errorbar(1:4,frac_card(1,:),sem(1,:),'k');hold on

end
end






    
