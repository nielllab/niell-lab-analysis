function barPinp(data,wt,ko1,ko2,lyr,pinped,inh,used)

% data_wt= nanmedian(data(wt & lyr&  ~pinped & ~inh & used));
% data_wt_sem=semedian(data(wt & lyr& ~pinped & ~inh & used));
[fWT,xWT] = ecdf(data((wt  & pinped  & used)));

data_WT=data((wt  & pinped  & used));
data_wt_p = nanmedian(data(wt  & pinped  & used));
data_wt_psem=semedian(data(wt  &  pinped  & used));

% data_ko= nanmedian(data(ko1 & lyr&  ~pinped & ~inh & used));
% data_ko_sem=semedian(data(ko1 & lyr& ~pinped & ~inh & used));
[fN2A,xN2A] = ecdf(data((ko1   & pinped  & used)))

data_N2A=data((ko1   & pinped  & used));
data_ko_p = nanmedian(data(ko1 & pinped  & used));
data_ko_psem=semedian(data(ko1  & pinped  & used));

% data_ko2= nanmedian(data(ko2 & lyr&  ~pinped & ~inh & used));
% data_ko2_sem=semedian(data(ko2 & lyr& ~pinped & ~inh & used));
[fN2B,xN2B] = ecdf(data((ko2   & pinped  & used)))

data_N2B=data((ko2  & pinped  & used))
data_ko_p2 = nanmedian(data(ko2  & pinped  & used));
data_ko2_psem=semedian(data(ko2  & pinped  & used));

% plot cumulative density plot
figure
plot(xWT,fWT,'g'); hold on
plot(xN2A,fN2A,'r'); hold on
plot(xN2B,fN2B,'b'); 

% plot bar graphs
figure
barweb([data_wt_p data_ko_p data_ko_p2],[ data_wt_psem data_ko_psem data_ko2_psem]);
%barweb([data_wt data_ko data_ko2 ; data_wt_p data_ko_p data_ko_p2],[data_wt_sem data_ko_sem data_ko2_sem ; data_wt_psem data_ko_psem data_ko2_psem]);
legend({'wt','N2A','N2B'})
set(gca,'xticklabel',{'pinp'});


% plot distributions in terms of fraction of total population
% figure 
% [f,x]=hist(data_WT,0.01:0.05:0.35);
% H1=bar(x,f/sum(f),'g');
% title 'SF pref WT'
% hold on
% 
% %figure
% [f1,x1]=hist(data_N2A,0.01:0.05:0.35);
% H2=bar(x1,f1/sum(f1),'r');
% title 'SF pref N2A'
% hold on
% 
% %figure
% [f2,x2]=hist(data_N2B,0.01:0.05:0.35);
% H3=bar(x2,f2/sum(f2),'b');
% title 'SF pref N2B'

 figure 
[f,x]=hist(data_WT);
H1=bar(x,f/sum(f),'g');
title 'SF pref WT'
hold on

%figure
[f1,x1]=hist(data_N2A);
H2=bar(x1,f1/sum(f1),'r');
title 'SF pref N2A'
hold on

%figure
[f2,x2]=hist(data_N2B);
H3=bar(x2,f2/sum(f2),'b');
title 'SF pref N2B'


figure
hist(data((wt  & pinped & ~inh & used)));
figure
hist(data((ko1   & pinped & ~inh & used)));
title 'N2A KO'
figure
hist(data((ko2   & pinped & ~inh & used)));
title 'N2B KO'
% 
% 







