function subplotAxes(gca,n)
set(gca,'Fontsize',10)
if ~exist('n','var') | n==5
   set(gca,'Xticklabel',{'sON','sOFF','tOFF','DS/OS','w-like'});
else
    
    set(gca,'Xticklabel',{'sON','sOFF','tOFF','DS/OS','w-like','supp'});
end
rotateticklabel(gca,90)
set(get(gca,'Ylabel'),'FontSize',10); 
