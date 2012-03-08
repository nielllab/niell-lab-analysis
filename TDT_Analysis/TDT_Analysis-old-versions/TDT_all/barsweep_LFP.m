% Matlab codes for reading from TTank
% and plot histgrams and rasters with each orientation
% Jianhua Cang 06-27-03

% connect to the tank and read the block


clear all;
pack
%close all
clear all
Tank_Name='052906'
Block_Name='revcheck2'
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

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end

invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger

Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time
 
CSDall=0;

figure
for orientation =1
%orientation =1;
clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
    Epocs=find(MyEpocs(1,:)==orientation); 
    trial_no = length(Epocs);
    Epocs_TS = MyEpocs(2,Epocs);
    
    invoke(TTX,'ResetFilters');
  i=orientation;
  for trial = 1:2
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, trial, 0); 
                                 % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
   end
    tranges = invoke(TTX,'GetValidTimeRangesV')
    eps = invoke(TTX,'GetEpocsExV','xTrg',0)    
    clear phi
    temp_filt = ones(1,20);
    
    for ch = 1:16
 
        allsweeps = invoke(TTX, 'ReadWavesOnTimeRangeV',  Event_Name_LFP,ch_map(ch));
        phi(ch,:) = mean(allsweeps(:,:),2);
        minimum(ch,:) = min(allsweeps(1:100,:));
        minimum2(ch,:) = min(allsweeps(200:300,:));
        all(ch,:,:) = allsweeps(1:100,:);
%         figure
%         plot(mean(allsweeps(1:200,1:10),2)');
%         hold on
%         plot(mean(allsweeps(1:200,50:60),2)','r');
    end

    freq = fft(phi');
   
    nyq = 1/(2*Sample_Interval)
    freqInt = nyq/(0.5*size(phi,2))
%     for i = 1:2
%         maxFreq =25*i;
%         maxfreqIndex = round(maxFreq/freqInt)
%         freq(maxfreqIndex-5:maxfreqIndex+8,:)=0;
%         freq(size(phi,2)-maxfreqIndex-8:size(phi,2)-maxfreqIndex+5,:)=0;
%     end
    maxFreq = 100;
    maxfreqIndex = round(maxFreq/freqInt);
    freq(maxfreqIndex:size(phi,2)-maxfreqIndex,:)=0;
% 
    figure
    plot(abs(freq(1:100,:)));
    phi2 = ifft(freq)';
%     figure
%     plot(real(phi2)');
%     figure
%     plot(imag(phi2)');
    phi = real(phi2);
%     phi(1,:) = phi(4,:);
%     phi(2,:)=phi(4,:);
%     phi(3,:)=phi(4,:);
    
    %phi = imfilter(phi,ones(1,10),'replicate','same');
    meanLFP = ones(16,1)*mean(phi);
    figure
    plot((phi(1:7,1:200)'));
    figure
    plot(phi(10:16,1:200)');
    
%     figure
%     plot((phi(:,1:200)-meanLFP(:,1:200))');
%     figure
%     plot((phi(:,600:800)-meanLFP(:,600:800))');
%     figure
%     plot((phi(9:16,1:200)-meanLFP(9:16,1:200))');
%     figure
%     plot((phi(9:16,700:900)-meanLFP(9:16,700:900))');
    figure
    plot((phi-meanLFP)');
%     figure
%     plot((phi(1:7,:)-meanLFP(1:7,:))');
    
%     CSD = (phi(1:10,:) + phi(7:16,:) - 2*phi(4:13,:))/(-9);    
%     figure
%     imagesc(CSD(:,1:150),[-0.5*10^-5 0.5*10^-5]);
%     colorbar
%     space_filter = [0.5 ; 1; 0.5]
%     space_filter = space_filter/sum(sum(space_filter));
%     CSD = imfilter(CSD,space_filter,'same');
%     figure
%     imagesc(CSD(:,1:300));
%     colorbar
    
       CSD = (phi(1:12,:) + phi(5:16,:) - 2*phi(3:14,:))/(-4);    
    newc = imresize(CSD,[60 size(CSD,2)],'bilinear');
    figure 
    imagesc(newc(:,1:150),[-0.25*10^-5 0.25*10^-5]);
%         figure 
%     imagesc(newc(:,750:950),[-0.25*10^-5 0.25*10^-5]);
%     space_filter = [0.5 ; 1; 0.5]
%     space_filter = space_filter/sum(sum(space_filter));
%     CSD = imfilter(CSD,space_filter,'same');
%     figure
%     imagesc(CSD(:,1:300));
%     colorbar 
    
    
    
        CSD = -1*(phi(1:14,:) + phi(3:16,:) - 2*phi(2:15,:));    
        figure
    imagesc(CSD(:,1:150),[-0.5*10^-5 0.5*10^-5]);
        
        newc = imresize(CSD,[70 size(CSD,2)],'bilinear');
    figure 
    imagesc(newc(:,1:150),[-0.75*10^-5 0.75*10^-5]);
    
%     figure
%     for i = 1:16;
%         subplot(8,2,i);
%         %figure
%         plot(phi(i,1:200)-meanLFP(i,1:200));
%         axis([0 200 -1*10^-5 1*10^-5])
%         axis off
%     end
    % figure
    % for i = 1:12;
    %     subplot(6,2,i);
    %     %figure
    %     plot(CSD(i,:));
    %     axis([700 900 -10^-3 10^-3])
    %     axis off
    % end
    
    %CSDall = CSDall+CSD;
end


    


invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');