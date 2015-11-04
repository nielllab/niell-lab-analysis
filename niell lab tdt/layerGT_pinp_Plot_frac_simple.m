function [meandata se N1 ]=layerGT_pinp_Plot_frac_simple(data,pinp,geno,layer,inh,exc,used,label);
colorlist='bmkg';
pinp=pinp;
geno=geno';

for GT=1:3
    for P=1:2
        
        for group = 1:1
            if group ==1
                uselist = (geno==GT & pinp==P & layer & used);
%             elseif group ==2
%                 uselist = (geno==GT & pinp==P & (layer==4) & exc & used);
%             elseif group==3
%                 uselist = (geno==GT & pinp==P & (layer==5) & exc & used);
            end
            
            lin = histc(data(uselist),0.8:1:2);
            countdata_lin(group,P) = lin(1,1);
            N(group,P) = sum(~isnan(data(uselist)));
            frac(group,P)= countdata_lin(group,P)/N(group,P);
            %       P_FF(group,age)= countdata_FF(group,age)/N(group,age)
            
            [fr,pci]= binofit(countdata_lin(group,P),N(group,P));
            errdata(group,P) = frac(group,P)-pci(1,1);
            sem(group,P)=errdata(group,P)/sqrt(countdata_lin(group,P));
            
        end
    end

% figure
% errorbar(1:3,frac(1,:),sem(1,:),'k');hold on
% errorbar(1:3,frac(2,:),sem(2,:),'g');hold on

% for group = 1:2
%  shadedErrorBar(1:4,frac(group,:),prct_err_lin(group,:),'k');
%     hold on;
% end
 
%  ylabel(label{1,1});
%  set(gca,'Xtick',1:4);
%  set(gca,'Xticklabel',{'EO1','EO3','EO7','Adult'});

% ranksum(bothdata{1,1},bothdata{1,2}) % signfiicance layer2/3
% ranksum(bothdata{2,1},bothdata{2,2}) % signfiicance layer4
% ranksum(bothdata{3,1},bothdata{3,2}) % sign. layer5
% ranksum(bothdata{4,1},bothdata{4,2}) % sign layer6
% % ranksum(bothdata{5,1},bothdata{5,2})% sign inhibotyr
% ranksum(bothdata{6,1},bothdata{6,2}) % all excitaor
%Chi_square(2,groups(1,1:2),groups(1,3:4));


figure
barweb(frac,sem)
ylabel(label);
set(gca,'Xtick',1:4);
set(gca,'Xticklabel',{'2/3','4','5'});
legend('Not-PINPed','PINPed');
title(GT)

end
end










    
