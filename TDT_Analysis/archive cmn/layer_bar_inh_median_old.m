function layer_bar_inh_median_old(data,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix)
figure
barweb([median(data(used&exc&layertwo)), median(data(used&exc&layerthree)), median(data(used&exc&layerfour)) ,median(data(used&exc&layerfive)),median(data(used&exc&layersix)),median(data(used&inh)); 0 0 0 0 0 0]', ...
    [std(data(used&exc&layertwo))/sqrt(sum(used&exc&layertwo )), std(data(used&exc&layerthree))/sqrt(sum(used&exc&layerthree)) ,  ...
    std(data(used&exc&layerfour))/sqrt(sum(used&exc&layerfour )), std(data(used&exc&layerfive))/sqrt(sum(used&exc&layerfive )) , std(data(used&exc&layersix))/sqrt(sum(used&exc&layersix)), std(data(used&inh))/sqrt(sum(used&inh)); 0 0 0 0 0 0]', ...
    2,{'upper 2/3' 'lower 2/3' '4' '5' '6' 'inh'},[],[],[],[0 0 .75]);
set(gca,'xticklabel', {'upper 2/3', 'lower 2/3', '4', '5' '6' 'inh'})
set(gca,'FontSize',14)

