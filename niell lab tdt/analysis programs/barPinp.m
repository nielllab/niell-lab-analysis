function barPinp(data,wt,ko1,ko2,lyr,pinped,inh,used)

% data_wt= nanmedian(data(wt & lyr&  ~pinped & ~inh & used));
% data_wt_sem=semedian(data(wt & lyr& ~pinped & ~inh & used));

data_wt_p = nanmedian(data(wt & lyr& pinped & ~inh & used));
data_wt_psem=semedian(data(wt & lyr& pinped & ~inh & used));

% data_ko= nanmedian(data(ko1 & lyr&  ~pinped & ~inh & used));
% data_ko_sem=semedian(data(ko1 & lyr& ~pinped & ~inh & used));

data_ko_p = nanmedian(data(ko1 & lyr& pinped & ~inh & used));
data_ko_psem=semedian(data(ko1 & lyr& pinped & ~inh & used));

% data_ko2= nanmedian(data(ko2 & lyr&  ~pinped & ~inh & used));
% data_ko2_sem=semedian(data(ko2 & lyr& ~pinped & ~inh & used));

data_ko_p2 = nanmedian(data(ko2 & lyr& pinped & ~inh & used));
data_ko2_psem=semedian(data(ko2 & lyr& pinped & ~inh & used));

figure
barweb([data_wt_p data_ko_p data_ko_p2],[ data_wt_psem data_ko_psem data_ko2_psem]);
%barweb([data_wt data_ko data_ko2 ; data_wt_p data_ko_p data_ko_p2],[data_wt_sem data_ko_sem data_ko2_sem ; data_wt_psem data_ko_psem data_ko2_psem]);
legend({'wt','N2A','N2B'})
set(gca,'xticklabel',{'pinp'});
