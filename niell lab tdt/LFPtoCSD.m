%function LFPtoCSD
% Code to compute CSD from LFP in response to 2-phase stimuli
% i.e. reversing checkerboard, grating, or full-field
% cmn 06-06

% connect to the tank and read the block
clear all;
pack
clear all;  % something in matlab memory management requires 2 clears ...???
  pname = uigetdir('C:\data\TDT','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
TTX = openTTX(Tank_Name,Block_Name); % to initialize
Event_Name_LFP = 'Lfpx'
Sample_Interval=10^-3 * 0.04096 *32 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; 

%%% remap site numbers if using a connector meant for tetrode configuration
%%% (not used since 9/05)
tetrode_linear=0;
if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:64;
end


%%% read in waveforms
invoke(TTX,'CreateEpocIndexing');
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
invoke(TTX,'ResetFilters');
for ph = 1:2    %%% read in both phases
    ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
    invoke(TTX,'SetFilter', ecode, 69, ph, 0); 
                             % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
end
invoke(TTX,'GetValidTimeRangesV')
invoke(TTX,'GetEpocsExV','xTrg',0)    
for ch = 1:16
    allsweeps = invoke(TTX, 'ReadWavesOnTimeRangeV',  Event_Name_LFP,ch_map(ch));
    phi(ch,:) = mean(allsweeps(:,:),2);
end


%%% temporal filtering of LFP (done in frequency domain)
maxFreq = 200;
freq = fft(phi');
nyq = 1/(2*Sample_Interval)
freqInt = nyq/(0.5*size(phi,2))
maxfreqIndex = round(maxFreq/freqInt);
freq(maxfreqIndex:size(phi,2)-maxfreqIndex,:)=0;
figure
plot(abs(freq(1:200,:)));
phi = real(ifft(freq)');

figure
plot((phi(1:16,1:200)'));

figure
plot((phi([1:14 16],1:100)'));

LFPfig=figure;
for c= 0:3   
    subplot(2,2,c+1);
    plot((phi((4*c+1):(4*c+4),1:200)'));
    axis([0 100 -20*10^-4 6*10^-4])
end
title(Tank_Name);

figure
subplot(1,3,1);
  plot(phi(3:6,1:100)');
    axis([0 100 -6.5*10^-5 6*10^-5])
 subplot(1,3,2);
  plot(phi(7:10,1:100)');
    axis([0 100 -6.5*10^-5 6*10^-5])
    subplot(1,3,3);
  plot(phi(11:14,1:100)');
    axis([0 100 -6.5*10^-5 6*10^-5])   
    
  
meanLFP = ones(16,1)*mean(phi);
figure
%plot((phi-meanLFP)');
plot(phi(:,1:150)');
title(Tank_Name);


%%% calculate CSD with 2-site spacing
CSD = (phi(1:12,:) + phi(5:16,:) - 2*phi(3:14,:))/(-4);
figure
imagesc(CSD(:,1:100),[-0.5*10^-3 0.5*10^-3]);   
newc = imresize(CSD,[600 5*size(CSD,2)],'bilinear');  %%% interpolate at 10x sampling
CSDfig =figure ;
imagesc(newc(:,1:500),[-5*10^-5 5*10^-5]);
title(Tank_Name);

% 
% %%% calculate CSD with 3-site spacing
% CSD = (phi(1:10,:) + phi(7:16,:) - 2*phi(4:13,:))/(-9);    
% newc = imresize(CSD,[60 size(CSD,2)],'bilinear');  %%% interpolate at 10x sampling
% figure 
% imagesc(newc(:,20:100),[-1*10^-5 1*10^-5]);

%%% same thing at 1-site spacing
CSD = -1*(phi(1:14,:) + phi(3:16,:) - 2*phi(2:15,:));    
figure
imagesc(CSD(:,1:100),[-0.5*10^-5 0.5*10^-5]);        
newc = imresize(CSD,[700 5*size(CSD,2)],'bilinear');
figure 
imagesc(newc(:,1:500),[-0.5*10^-5 0.5*10^-5]);
title(Tank_Name);

resp=input('save results y/n','s')
if resp=='y'
    output_path=uigetdir('','data folder');
%          fname = fullfile(output_path,sprintf('CSD%s_%s',Tank_Name,Block_Name));
%         saveas(CSDfig,fname,'fig');
        fname = fullfile(output_path,sprintf('tetLFP%s_%s',Tank_Name,Block_Name));
        saveas(LFPfig,fname,'fig');

end

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');