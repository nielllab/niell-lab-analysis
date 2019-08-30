function corrPlot(a,b,c,d,labels,figTitle)
%%% make correlation figure plots to go in compileCap
%%% call as e.g. 
%%% labels {'left theta','right theta','left phi','right phi'}
%%% corrPlot(corrRAll,corrLAll,corrRAAll,corrLAll,labels,'all eye correlations')

data{1} = a; data{2} = b; data{3} =c; data{4}=d;
figure('units','normalized','outerposition',[0 0 1 1]); hold on

col = {'b-','r-','g-','c-'};
for i = 1:4
    err{i} = nanstd(data{i})/sqrt(length(data{i}));
    shadedErrorBar(1:size(data{i},2),nanmean(data{i},1),err{i},col{i},1);
end

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([0 1]); 
xlim([21 41]); axis square
for i = 1:4
    L(i) = plot(nan, nan, col{i});
end
legend(L,labels); title(figTitle);
