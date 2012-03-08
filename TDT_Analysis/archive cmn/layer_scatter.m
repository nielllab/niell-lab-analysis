function layer_scatter(data,used,layer, exc, inh);
figure
plot(data(used&inh),layer(used&inh)+rand(sum(used&inh),1)*0.4-0.2,'bo','MarkerSize',8,'LineWidth',1);
hold on
plot(data(used&exc),layer(used&exc)+rand(sum(used&exc),1)*0.4-0.2,'gs','MarkerSize',8,'LineWidth',1);
axis([min(data(used)) max(data(used)) 1.5 6.5])
axis ij
legend('narrow','broad')
set(gca,'ytick',[2 3 4 5 6])
set(gca,'FontSize',14)
ylabel('layer','FontSize',16)