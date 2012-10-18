afile = { ...%'C:\data\ephys matlab data\022012_awake_pptg\mlr1\analysis.mat',...
    'C:\data\ephys matlab data\022012_awake_pptg\mlr_rec2\analysis.mat','C:\data\ephys matlab data\070512_awake_mlr\MLRrec2a\analysis.mat',...
    'C:\data\ephys matlab data\070512_awake_mlr\MLRrec3\analysis.mat','C:\data\ephys matlab data\070712_awake_mlr\MLRrec1\analysis.mat',...
    'C:\data\ephys matlab data\071212_awake_mlr\MLRrec5_nostim\analysis.mat','C:\data\ephys matlab data\071212_awake_mlr\MLRrec6_nostim\analysis.mat',...
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
h = hist([co_all ; shuffle_co_all]',-0.5:0.25:0.75)
bar( -0.5:0.25:0.75, h/length(co_all))
legend({'raw','shuffled'})
xlabel('correlation')
xlim([-1 1])
ylabel('fraction')

figure
plot(d_all');
 