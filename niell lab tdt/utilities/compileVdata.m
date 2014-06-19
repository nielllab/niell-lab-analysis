clear all
afile = {...
%    'C:\data\ephys matlab data\021412_awake_chr2\wndrift1\analysis.mat', ...     %%LFP problems
%     'C:\data\ephys matlab data\021512_awake_pptg\wn1b\analysis.mat', ...        %%LFP problems
     'C:\data\ephys matlab data\021612_awake\wn1drift1\analysis.mat', ...         % MLRstim #1; great desynchronization
    'C:\data\ephys matlab data\021612_awake\wn2b\analysis.mat',...                % MLRstim #2; good desychronization
    'C:\data\ephys matlab data\021612_awake\wn3\analysis.mat',...                 % MLRstim #3; good desychronization at 50Hz
    'C:\data\ephys matlab data\021612_awake\wn4bdrift5\analysis.mat',...          % MLRstim #4    
    'C:\data\ephys matlab data\021812_awake_pptg\wn2\analysis.mat',...            % MLRstim #5    
    'C:\data\ephys matlab data\021812_awake_pptg\wn3e_drift2\analysis.mat',...    % MLRstim #6 % bad lfp
    'C:\data\ephys matlab data\021812_awake_pptg\wn5\analysis.mat',...            % MLRstim #7
    'C:\data\ephys matlab data\070112_awake_mlr\wn1fg\analysis.mat',...           % MLRstim #8  % bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn2a\analysis.mat',...            % MLRstim #9  % bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn3d\analysis.mat',...            % MLRstim #10
    'C:\data\ephys matlab data\070112_awake_mlr\wn4\analysis.mat',...             % MLRstim #11  %%% bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn5b\analysis.mat',...            % MLRstim #12
    'C:\data\ephys matlab data\070312_awake_mlr\wn1h\analysis.mat',...            % MLRstim #13
    'C:\data\ephys matlab data\070312_awake_mlr\wn2e\analysis.mat',...            % MLRstim #14
    'C:\data\ephys matlab data\070312_awake_mlr\wn3a\analysis.mat',...            % MLRstim #15
   'C:\data\ephys matlab data\070312_awake_mlr\wn4c\analysis.mat'...              % MLRstim #16
}

% afile = {'C:\data\ephys matlab data\070612_awake_mlr\wn4b\analysis.mat',...    %MLR2BFstim#1
%     'C:\data\ephys matlab data\070612_awake_mlr\wn5b\analysis.mat',...      %MLR2BFstim#2
%     'C:\data\ephys matlab data\070612_awake_mlr\wn6b\analysis.mat',...      %MLR2BFstim#3
%     'C:\data\ephys matlab data\070612_awake_mlr\wn7a\analysis.mat',...      %MLR2BFstim#4
%     'C:\data\ephys matlab data\070612_awake_mlr\wn8b\analysis.mat',...      %MLR2BFstim#4
%     'C:\data\ephys matlab data\070612_awake_mlr\wn9a\analysis.mat',...      %MLR2BFstim#5
%     'C:\data\tdt tanks\071312_awake_MLR\wn1b\analysis.mat',...              %MLR2BFstim#6
%     'C:\data\ephys matlab data\071312_awake_mlr\wn2d\analysis.mat',...      %MLR2BFstim#7
%     'C:\data\ephys matlab data\071312_awake_mlr\wn5\analysis.mat',...       %MLR2BFstim#8
%     'D:\Moses_ephys_data\082012_awake_mlr\wn2a\analysis.mat',...            %MLR2BFstim#9
%     'D:\Moses_ephys_data\082012_awake_mlr\wn3b\analysis.mat',...            %MLR2BFstim#10
%     'D:\Moses_ephys_data\082112_awake_mlr\wn2c\analysis.mat',...            %MLR2BFstim#11
%     'D:\Moses_ephys_data\082312_awake_MLR\wn3d\analysis.mat',...            %MLR2BFstim#12
%     'D:\Moses_ephys_data\082312_awake_MLR\wn5a\analysis.mat',...            %MLR2BFstim#13
%     'D:\Moses_ephys_data\082312_awake_MLR\wn6b\analysis.mat'}               %MLR2BFstim#14
%     

vdataAll=[];
vMeanAll=[];
vFracAll=[];
for i = 1:length(afile)
    load(afile{i});
    vdataAll = [vdataAll vdata];
    vMeanAll = [vMeanAll mean(vdata,2)];
    vFracAll = [vFracAll mean(vdata>1,2)];
end


figure
plot(0:0.1:35,vFracAll,'Color',[0.7 0.7 0.7]);
hold on
plot(0:0.1:35,mean(vFracAll,2),'b','Linewidth',4);
plot([10 27],[0.9 0.9],'g','Linewidth',8)
xlabel('secs');
ylabel('% running');



figure
plot(0:0.1:35,median(vdataAll,2));
hold on
plot(0:0.1:35,prctile(vdataAll,75,2));
plot(0:0.1:35,prctile(vdataAll,25,2));
xlabel('secs');
ylabel('cm/sec');
ylim([0 1]);
hold on
plot([10 27],[4.9 4.9],'g','Linewidth',8)


figure
plot(0:0.1:35,mean(vdataAll>1,2));
xlabel('secs');
ylabel('% time running');
ylim([0 1]);
hold on
plot([10 27],[0.9 0.9],'g','Linewidth',8)
title('NB term stim')


nonstim = [1:100 271:350];
stim = 101:270;
fracnon = mean(vFracAll(nonstim,:),1)
fracstim = mean(vFracAll(stim,:),1)

figure
title('NB term stim')
barwitherr([std(fracnon)/sqrt(length(fracnon)) std(fracstim)/sqrt(length(fracstim))], ...
    [mean(fracnon) mean(fracstim)]);
ylabel('% time running');
ylim([0 1]);
set(gca,'XTickLabel',{'pre/post stim','stim'})
