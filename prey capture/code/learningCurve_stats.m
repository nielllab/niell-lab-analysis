
MED=nanmedian(timeToC)
sMed=semedian(timeToC)

figure
errorbar(1:5,MED,sMed,'k');hold on

figure
errorbar(1:3,MED(1,3:5),sMed(1,3:5),'k');hold on


clear n p
n=length(successCap);
p=nanmean(successCap);
q=1-p;

[m v]=binostat(n,p);
m=m/n; v=v/n;
stdDev_suc=sqrt(v);
SEM=stdDev_suc/sqrt(n);

figure
errorbar(1:5,m,SEM,'k');hold on


[tbl, chi2stat,pval]=crosstab(successCap(:,1),successCap(:,3));

%% plexiglass
clear n p
n=length(successCap);
p=nanmean(successCap);
q=1-p;

[m v]=binostat(n,p);
m=m/n; v=v/n;
stdDev_suc=sqrt(v);
SEM=stdDev_suc/sqrt(n);

figure
errorbar(1:5,m,SEM,'k');hold on


[tbl, chi2stat,pval]=crosstab(successCap(:,1),successCap(:,3));