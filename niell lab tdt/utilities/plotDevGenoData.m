function plotDevData(vals, age, agelist, ageBins, used, geno,label)

col = 'gbk';
figure
for g=0:2
[resp resperr] = sortbyage(vals,age,ageBins,used& geno==g)

hold on
errorbar(mean(ageBins(:,~isnan(resp)),1) + (g-1)*0.2,resp(~isnan(resp)), resperr(~isnan(resp)),[col(g+1) '-o']);
%plot(mean(ageBins(:,~isnan(resp)),1)+ (g-1)*0.2,resp(~isnan(resp)),[col(g+1) 'o']);
end

ylabel(label)
legend('ko','het','wt')