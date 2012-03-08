% Matlab codes for reading from TTank
% and plot histgrams and rasters with each orientation
% Jianhua Cang 06-27-03

% connect to the tank and read the block


clear all;
pack
%close all
clear all
Tank_Name='052308_awake_tet'
Block_Name='bars16d2'
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Event_Name_LFP = 'Lfpx'
Sample_Interval=10^-3 * 0.04096 *32 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=8; %in second
hist_range=[0:0.1:9];
axis_range=[0 plot_duration 0 4];


tetrode_linear=0;


invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger


 
CSDall=0;

nrows=5;
ncols=4;
ncond = 17;
ch=15;
wv = figure;
spec = figure;
figure
for orientation =1:ncond
%orientation =1;
clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
    Epocs=find(MyEpocs(1,:)==orientation); 
    trial_no = length(Epocs);
    Epocs_TS = MyEpocs(2,Epocs);
    
    invoke(TTX,'ResetFilters');
  i=orientation;

        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, orientation, 0); 
                                 % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'

    tranges = invoke(TTX,'GetValidTimeRangesV')
    eps = invoke(TTX,'GetEpocsExV','xTrg',0)    
    clear phi
    temp_filt = ones(1,20);

        allsweeps = invoke(TTX, 'ReadWavesOnTimeRangeV',  Event_Name_LFP,ch);
       figure(wv)
       subplot(nrows,ncols,orientation);
        imagesc(allsweeps',[-5*10^-4 5*10^-4]);
    figure(spec)
        subplot(nrows,ncols,orientation);
        onset_time =0.02;
        freq = size(allsweeps,1)/ (tranges(2,1)-tranges(1,1))
       allsweeps_crop = allsweeps(onset_time*freq:size(allsweeps,1),:);
        pwr = abs(fft(allsweeps_crop,[],1))./size(allsweeps_crop,1);
       compression_factor =10;
       pwr_short = condenseData(pwr,10);
       freq = size(allsweeps,1)/ (tranges(2,1)-tranges(1,1))
       allsweeps_crop = allsweeps(onset_time*freq:size(allsweeps,1),:);
       df = compression_factor / (tranges(2,1)-tranges(1,1)-onset_time);
       for i = 1:size(pwr_short,2);
           pwr_short(:,i)=pwr_short(:,i).*(1:size(pwr_short,1))'*df;
       end
       
        imagesc(pwr_short(1:(100/df),:)',[0 3*10^-4]);
       axis off
         set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
end




invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');