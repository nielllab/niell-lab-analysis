function scatterDeveGenoData(vals, age, agelist, ageBins, used, geno,label)

col = 'gbk';
figure
hold on
for g=0:2
plot(age(used&geno==g)+ 0.3*(g-1) + rand(size(age(used&geno==g)))*.2,vals(used&geno==g),[col(g+1) 'o'])
[col(g+1) 'o']
end

ylabel(label)
legend('ko','het','wt')