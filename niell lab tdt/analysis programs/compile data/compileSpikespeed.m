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
end
    
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
 