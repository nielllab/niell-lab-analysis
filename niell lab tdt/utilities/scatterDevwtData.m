function scatterDevwtData(vals, age, agelist, ageBins, used, geno,label,splot)

col = 'gbk';

if ~exist('splot','var') | ~splot
    figure
end

hold on
%age(age>25)=25;
for g=0:2
plot(age(used&geno==g) + rand(size(age(used&geno==g)))*.2,vals(used&geno==g),[col(g+1) 'o'])
[col(g+1) 'o']
end

ylabel(label)

% set(gca,'Xtick',[14 16 18 20 22 28]);
% set(gca,'Xticklabel',{'14','16','18','20','22','adult'});
% xlabel('age');

set(gca,'Xtick',[14 16 18 20 22 24 28]);
set(gca,'Xticklabel',{'14','16','18','20','22','24','adult'});
xlabel('age');
xlim([13.5 29])