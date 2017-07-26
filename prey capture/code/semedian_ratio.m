function s= semedian_ratio(data,data1)
x=[data, data1];

s = nanstd(bootstrp(1000,@(x) ((nanmedian(x(:,1))-nanmedian(x(:,2)))/(nanmedian(x(:,1))+ nanmedian(x(:,2)))),x)) ;

%s = nanstd(bootstrp(1000,@(x) nanmedian(x,1),data));
% n=25;                   %number of neurons
% nReps = 10000;          %number of iterations for the bootstrap
% CIrange = 95;           %confidence interval range
% 
% x = ceil(15*randn(n,2).^2);      %nx2 matrix of firing rates (Chi-squared distribution)
% define our 'index' here (difference over sum)
% 
% myStatistic = @(x) mean((x(:,1)-x(:,2))./(x(:,1)+x(:,2)));
% run the 'boostrap' program to generate the confidence interval
% 
% [CI,sampStat,bootstrapStat] = bootstrap(myStatistic,x,nReps,CIrange);