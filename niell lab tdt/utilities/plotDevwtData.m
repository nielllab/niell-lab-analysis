function plotDevwtData(vals, age, agelist, ageBins, used, geno,label)

col = 'gbk';
figure

g=2
[resp resperr] = sortbyage(vals,age,ageBins,used& geno==g)

hold on
errorbar(mean(ageBins(:,~isnan(resp)),1) , resp(~isnan(resp)),resperr(~isnan(resp)),[col(g+1) '-o']);
%plot(mean(ageBins(:,~isnan(resp)),1)+ (g-1)*0.2,resp(~isnan(resp)),[col(g+1) 'o']);


ylabel(label)
set(gca,'Xtick',[14 16 18 20 22 24 28);
set(gca,'Xticklabel',{'14','16','18','20','22','24','adult'});
xlabel('age');
xlim([13.5 29])
