figure
for c=1:length(all_fit)
subplot(4,5,c)
imagesc(all_fit{c});axis square
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
end