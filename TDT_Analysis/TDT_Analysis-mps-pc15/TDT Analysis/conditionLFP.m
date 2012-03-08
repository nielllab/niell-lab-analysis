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
    Tank_Name = '032207_wt_tet';
    Block_Name = 'bars16d1c';
    time1= 0;
    time2=0;
    max_events = 50000;
end




TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_LFP='Lfpx'

%%% set time based on first and last timepoints in clustered data
 for cond = 16:17
    invoke(TTX,'CreateEpocIndexing');
    MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000)
    invoke(TTX,'ResetFilters');

    %%% read in both phases
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, cond, 0); 
                                 % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'

    invoke(TTX,'GetValidTimeRangesV')
    invoke(TTX,'GetEpocsExV','xTrg',0)    

    for ch = 1:16
        allsweeps = invoke(TTX, 'ReadWavesOnTimeRangeV',  Event_Name_LFP,ch);
        all_fft = fft(allsweeps(:,:),[],1);
        all_fft(1:90,:)=0;
        all_fft(161:1060,:)=0;
        all_fft(1130:size(all_fft,1),:)=0;
        allsweeps_filt = ifft(all_fft,[],1);
%        
%         if ch==9
%         figure
%         for i = 1:min(16,size(allsweeps_filt,2));
%               subplot(4,4,i);
%               %figure
%               plot(real((allsweeps_filt(:,i))));
%         end
%         end
        phi(ch,:) = mean(abs(fft(allsweeps(:,:),[],1)),2);

       phi_mean(ch,:) = mean(allsweeps(:,:),2);
        % phi_norm(ch,:) = phi(ch,:).*(1:size(phi,2));
    end


    figure
    phi_condense=condenseData(phi',5)';
for c = 1:16
    phi_condense(c,:) = phi_condense(c,:).*(1:size(phi_condense,2));
end
    plot(phi_condense(:,1:50)');
    axis([0 50 0 0.2])
%     figure
%     plot(phi_mean');

 end

% 
% %%% reshape wave data so that one FFT block runs down a column
% blk_size = 1900; %%% approximately 10 sec
% nblk = floor(size(W,1)/blk_size)
% 
% W_trim = W(1:blk_size*nblk);
% W_trim = reshape(W_trim,blk_size,nblk);
% 
% %%% take FFT and diplay it
% EEG_power = fft(W_trim,[],1);
% compression_factor = 10
% EEG_condense = condenseData(abs(EEG_power),compression_factor);
% freq = 32/(Wave_TS(5)-Wave_TS(4));
% freq_int = compression_factor*freq/blk_size   %%% this gives frequency scale for x-axis, not sure how to relabel the figure in imagesc
% figure
% imagesc(abs(EEG_condense(1:(100/freq_int),:))',[0 10*10^-3]);
% 
% spectrum = mean(EEG_condense(1:(100/freq_int),:),2);
% figure
% plot(spectrum);
% saveas(gcf,fullfile(pname,sprintf('EEG%s',Block_Name)),'fig');

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');