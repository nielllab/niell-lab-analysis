%function analyzeEEG
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block
clear all
clusterfile=0; %% use cluster file to get info, or not

if clusterfile
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    times = event_times_all(:);
    times = times(times>0);
    time1 = min(times);
    time2= max(times);
    block = input('Block number : ')    
    Block_Name = char(Block_Name(block));
else
    Tank_Name = '020109_awake_tet';
    Block_Name = 'bars16d1';
    time1= 0;
    time2=0;
    max_events = 50000;
end




TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_EEG='Lfpx'

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
nEpocs = size(MyEpocs,2);
 
%%% read in EEG wave info
% N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG, 10, 0, ...
%             MyEpocs(2,1),MyEpocs(2,nEpocs),'All');
N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,5, 0, ...
         0,0,'All');
        
if (N==max_events)
    warning('max number of events aquired');
end
W = invoke(TTX, 'ParseEvV', 0, N);
Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
W = W(:);

%%% reshape wave data so that one FFT block runs down a column
%%% samples at 381 Hz;
onesec=381;
duration = 4;
for i = 1:floor(length(W)/onesec)-duration;
    window(i,:) = W((i-1)*onesec+1:(i+duration-1)*onesec);
end
    
blk_size = duration*onesec; %%% approximately 5 sec
nblk = floor(size(W,1)/blk_size)

W_trim = W(1:blk_size*nblk);
W_trim = reshape(W_trim,blk_size,nblk);

W_trim = window';
%%% take FFT and diplay it
EEG_power = fft(W_trim,[],1);
compression_factor = 2
EEG_condense = condenseData(abs(EEG_power),compression_factor);
freq = 32/(Wave_TS(5)-Wave_TS(4));
freq_int = compression_factor*freq/blk_size   %%% this gives frequency scale for x-axis, not sure how to relabel the figure in imagesc
% figure
% %imagesc(abs(EEG_condense(1:(180/freq_int),:))',[0 2*10*10^-3]);
% imagesc(abs(EEG_condense(1:(180/freq_int),:))');

EEG = abs(EEG_condense(1:(75/freq_int),:))';
EEG_norm = zeros(size(EEG));
size(EEG_norm)
for i = 1:size(EEG,1);
    EEG_norm(i,:) = EEG(i,:).*(1:size(EEG,2));
end
figure
%imagesc(EEG_norm,[0 0.1]);
imagesc(EEG_norm);
set(gca,'XTick',(10:10:60)/freq_int);
set(gca,'XTickLabel',{'10','20','30','40','50','60'})


spectrum = mean(EEG_condense(1:(70/freq_int),:),2);
figure
plot(spectrum.*(1:size(spectrum,1))','r');
set(gca,'XTick',(10:10:60)/freq_int);
set(gca,'XTickLabel',{'10','20','30','40','50','60'})

figure
plot(mean(EEG_norm(:,ceil(7/freq_int):ceil(10/freq_int)),2));
hold on
plot(mean(EEG_norm(:,ceil(12/freq_int):ceil(18/freq_int)),2)+1,'g');
plot(mean(EEG_norm(:,ceil(20/freq_int):ceil(30/freq_int)),2)+2,'r');
plot(mean(EEG_norm(:,ceil(55/freq_int):ceil(62/freq_int)),2)+3,'k');



saveas(gcf,fullfile(pname,sprintf('EEG%s',Block_Name)),'fig');

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');