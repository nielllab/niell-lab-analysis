clear all
close all
afile = {...%'C:\data\ephys matlab data\022012_awake_pptg\mlr1\analysis.mat',...
    'C:\data\ephys matlab data\022012_awake_pptg\mlr_rec2\analysis.mat', ...   %updated
    'C:\data\ephys matlab data\070512_awake_mlr\MLRrec2a\analysis.mat',... %%% updated
    'C:\data\ephys matlab data\070512_awake_mlr\MLRrec3\analysis.mat', ... %%%updated
    'C:\data\ephys matlab data\070712_awake_mlr\MLRrec1\analysis.mat',... %%% updated
    'C:\data\ephys matlab data\071212_awake_mlr\MLRrec5_nostim\analysis.mat', ...  %%%updated
    'C:\data\ephys matlab data\071212_awake_mlr\MLRrec6_nostim\analysis.mat',... %%%updated
    'D:\Moses_ephys_data\082112_awake_mlr\MLRrec1b\analysis.mat'};

hist_all=0;
hist_all_shuffle=0;
n=0;
for i= 1:length(afile);
    load(afile{i})
        cell_range = n+(1:length(co));
    n=max(cell_range);
    
    co_all(cell_range) = co;
    shuffle_co_all(cell_range) =shuffle_co;
    
    d_all(cell_range,:)=d;
    for j=1:length(cell_range)
        R{cell_range(j)} = Rall{j};
    mousehist{cell_range(j)} = mouseVall{j};
    end
end
   


for cell_n = 1:n

    
    [c,pvals] = corrcoef(mousehist{cell_n},R{cell_n});
    [rankc ] = corr(mousehist{cell_n}',R{cell_n}','type','Spearman')
    % title(sprintf('corr = %f rho =%f',c(1,2),rankc))
     co(cell_n) = c(1,2);
     allPvals(cell_n)=pvals(1,2);
     
     shifts=-2:2;
     xc=zeros(size(shifts));
     for shift=1:length(shifts);
         xc(shift)=corr(circshift(mousehist{cell_n}',shifts(shift)),R{cell_n}','type','Spearman');
     end
     figure
     plot(shifts*2,xc); xlabel('lag (secs)'); ylabel('rank corr');
     
     [m ind] = max(abs(xc));
     maxC(cell_n) = xc(ind);
     maxT(cell_n) = shifts(ind)*2;
     title(sprintf('maxC = %f maxT = %f',maxC(cell_n),maxT(cell_n)));
     
     shiftC=zeros(length(R)-1,1);
     for shift=1:length(R)-1;
         shiftC(shift) = corr(circshift(mousehist{cell_n}',shift),R{cell_n}','type','Spearman');
     end
%      figure
%      hist(shiftC)
     %%% calculate prctile and p-value
     
     
     
     %rankCorr(cell_n) = rankc;
      rankCorr(cell_n) = maxC(cell_n);
     if rankCorr(cell_n)<prctile(shiftC,2.5) | rankCorr(cell_n)>prctile(shiftC,97.5)
    sigCorr(cell_n) = 1;
     else
         sigCorr(cell_n)=0;
     end
     rankpctile(cell_n) = length(find(shiftC<rankCorr(cell_n)))/length(shiftC);
     
end


h_range=-1:0.2:1;
figure
h=hist(rankCorr,h_range);
bar(h_range,h/length(rankCorr),'w');
hold on
h=hist(rankCorr(find(sigCorr)),h_range);
bar(h_range,h/length(rankCorr));
xlim([-1 1]);
ylim([0 0.3])
xlabel('correlation');
ylabel('fraction');


figure
pie([ sum(rankCorr>0 & sigCorr==1) sum(rankCorr<0 & sigCorr==1)  sum(sigCorr==0)],{'locomotor-active','locomotor-suppressed','not significant'})

figure
bar([ sum(rankCorr>0 & sigCorr==1) sum(rankCorr<0 & sigCorr==1)  sum(sigCorr==0)]/length(sigCorr))
set(gca,'Xticklabel',{'locomotor-active','locomotor-suppressed','not significant'});
ylabel('fraction');
keyboard

figure
%%
h1 = hist(co_all',-1:0.1:1);
h2 = hist(shuffle_co_all',-1:0.1:1);

bar( -1:0.1:1, h1/length(co_all),'b'); hold on;
bar( -1:0.1:1, h2/length(co_all),'r');

legend({'raw','shuffled'})
xlabel('correlation')
xlim([-.8 .9])
ylabel('fraction')

loco_active=sum(co_all>.2)/length(co_all)
loco_suppress=sum(co_all<-.2)/length(co_all)
unclassified=1-loco_active-loco_suppress

figure
pie([loco_active loco_suppress unclassified],{'locomotor-active','locomotor-suppressed','unclassified'});

% figure
% plot(d_all');
 