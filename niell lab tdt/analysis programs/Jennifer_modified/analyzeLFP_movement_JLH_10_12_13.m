%function analyzeLFP_movement
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block
% 
% [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
% %psfilename = 'c:/test.ps';   %%% default location
% if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


% pname = uigetdir('C:\data\tdt tanks','block data')
close all

pname = 'D:\Jen tanks\7_17_13\mouseA_bars_1';
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan = 12; % input('number of chans : ');
movement = 1; % input('movement data 0/1 : ');

tic
if movement==0
    flags =struct('lfpTseries',1,'lfpSpectra',1);
else
    flags =struct('lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
end

keyboard

tdtData= getTDTdata(Tank_Name, Block_Name, nChan, flags);

keyboard

for ch = nChan;
   
    %load and normalize LFP data
    lfp = double(tdtData.lfpData{ch}); % tdtData.spectData{ch}; %cast is necessary to run filter and so that chronux gives us the same frequency buckets
    
        Fs = 1/median(diff(tdtData.lfpT{ch}));

%     lfp = randn(size(lfp)) + sin(2*pi*tdtData.lfpT{ch}'*60);
    
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
    
%keyboard    
    Wo = 60/(Fs/2);
    BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    
%     b = fir1(300,Wo+BW*[-1 1]/2,'stop');
%     a = 1;
    
    lfp_filter = filtfilt(b,a,lfpnorm);
    
    
    %[lfpT, lfpData, spectT, specF ,lfpnorm   ] = analyzeLFP_chronux([tdtData.lfpT{ch} lfpnorm'   ],ch,true,true);
    [lfpT, lfpData, spectT, specF2,lfp_filter] = analyzeLFP_chronux([tdtData.lfpT{ch} lfp_filter'],ch,true,true);

    if ~all(specF{ch}==specF2{ch})
        error('huh?')
    else 
        specF=specF{ch};
    end
    
    lfpnorm = lfpnorm{ch};
    lfp_filter = lfp_filter{ch};
%    keyboard
    
%     
%     lfp_filter = abs (fft (lfp_filter));
%     lfpnorm = abs (fft (lfpnorm));
    
%     fs = 762.94;
%     F0 = 60;
%     fn = fs/2;
%     freqRatio = F0/fn;
%     
%     notchWidth = .5
%     
%     zeros = [exp( sqrt(-1)*pi*freqRatio),exp(-sqrt(-1)*pi*freqRatio)];
%     
%     poles = (1-notchWidth)*zeros
%     
%     figure;
%     zplane(zeros.',poles.');
%     b = poly(zeros);
%     a = poly(poles);
%     
%     
%     figure;
%     freqz(b,a,32000,fs)
%     
%     lfp_filter = filter(b,a,lfpnorm);
    
 
    
%     %filter 60Hz
%     N     = 6;       % Order
%     F0    = 60;      % Center frequency
%     BW    = 6;       % Bandwidth
%     Astop = 100;     % Stopband Attenuation (dB)
%     Fs    = 762.94;  % Sampling Frequency
%     
%     h = fdesign.notch('N,F0,BW, N, F0, BW, Astop, Fs, 'Squared');
%     
%     Hd = design(h, 'butter', ...
%         'SOSScaleNorm', 'Linf');
%     
%     set(Hd, 'Arithmetic', 'double');
%     
%     lfp_filter = filter(Hd,lfpnorm);


%     H = fspecial('average',[3 3]) %smoothing interpolated data
%     lfp_filter = imfilter(lfp_filter,H);
%      
    
    
    figure
    imagesc(lfp_filter',[0 prctile(lfp_filter(:),95)]);
    axis xy
    df = median(diff(tdtData.spectF{ch}));
    dt = median(diff(tdtData.spectT{ch}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    
%        
%     if movement
%         hold on
%         tsamp = tdtData.mouseT;
%         vsmooth = tdtData.mouseV;
%         
%         %plot(tsamp,(vsmooth/1.3-40),'g');
%         plot(tsamp,(vsmooth/.2-40),'g');
%         axis([0 max(tsamp) -40 80/df])
%         
%         set(gcf, 'PaperPositionMode', 'auto');
%         print('-dpsc',psfilename,'-append');
%     end
%     title(sprintf('channel = %d',ch));
%     
   
    %%%%

    theta = mean(lfpnorm(:,ceil(7/df):ceil(10/df)),2);
    gamma = mean(lfpnorm(:,ceil(50/df):ceil(58/df)),2);
    
    %v_interp = interp1(tsamp,vsmooth,tdtData.spectT{ch});

    %     figure
    %     plot(v_interp,gamma(t),'o');
    %     figure
    %     plot(v_interp,theta(t),'o');
    Smean = mean(lfp_filter,2)';
%     stationary = find(v_interp<0.3 & Smean<(5*median(Smean)));
%     moving = find(v_interp>0.35  & Smean<(5*median(Smean)));

    figure
    plot(specF,mean(lfp_filter,1));
    hold on
    plot(specF,mean(lfpnorm,1),'g');
%     plot(mean(lfp_filter(moving,:),1),'g');
if false
    axis([0 70/df 0 1.2*max(mean(lfp_filter))]);
    set(gca,'XTick',(10:10:70)/df);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','70'})
end    
%     title(sprintf('site %d',ch));
%      set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');

    %close all
end %% tet
% 
% ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
% delete(psfilename);


