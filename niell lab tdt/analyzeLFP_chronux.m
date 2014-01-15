function [lfpT lfpData spectT spectF spectData] = analyzeLFP_chronux(TTX,chans,timeSeriesOn,spectrumOn);
% Matlab code for reading LFP data and
% performing running spectrum analysis with chronux
% cmn 06/2011

%oldFormat = input('old TDT format (pre-UO) ? 0/1 : ')
%%% need to switch this to read in old UCSF data
oldFormat = 0;
if oldFormat == [];
    oldFormat = 0;
end
if oldFormat
    Event_Name='Lfpx'
else
    Event_Name='pLFP'
end
max_events = 10^6;
max_t = 10^9;
LFPall=0;
for ch = chans
    
    [lfp ts] = readWave(TTX,ch, Event_Name,max_events,max_t);
    lfpd = double(lfp); %modified to be able to filter 60Hz, double is applied so that filtfilt may be applied to proper data type
    %filter 60Hz noise
    
    notchFilter=0;   %%% select
    if notchFilter
        Fs = 1/median(diff(ts));%762.936;
        
        Wo = 60/(Fs/2); %JLH
        BW = Wo/35; %JLH
        [b,a] = iirnotch(Wo,BW); %JLH
        
        lfp_filter = filtfilt(b,a,lfpd); %JLH
        
    else
        lfp_filter = lfp;
    end
    
    
    if timeSeriesOn
        lfpT{ch} = ts;
        lfpData{ch} = lfp_filter; %lfp_filter was lfp
        
     
    end
    if spectrumOn
        
        
        
        params.Fs = 1/median(diff(ts));
        params.tapers = [3 5];
        params.fpass = [0 150];
        
        % lfp=rmlinesmovingwinc(lfp,[1 1],[],params,[],'y');
        
        % lfp=rmlinesc(lfp,params,[],'y',60);
        
        
        [S t f] = mtspecgramc(lfp_filter',[3 1],params); %lfp_filter was lfp
        %[S t f] = mtspecgramc(lfp',[10 2],params);
        spectData{ch} = S;
        spectT{ch} = t;
        spectF{ch} = f;
    end
    %     figure
    %     imagesc(Snorm',[0 10^-7]);
    %     axis xy
    %     set(gca,'YTick',(10:10:80)/df);
    %     set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    
    %
    %         figure
    %         %imagesc(Snorm,[0 0.1]);
    %         imagesc(Snorm',[0 6*10^-8]);
    %         axis xy
    %         set(gca,'YTick',(10:10:80)/df);
    %         set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    
end %% tet