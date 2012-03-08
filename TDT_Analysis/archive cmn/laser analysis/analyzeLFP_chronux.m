function [lfpT lfpData spectT spectF spectData] = analyzeLFP_chronux(TTX,nChan,timeSeriesOn,spectrumOn);
% Matlab code for reading LFP data and 
% performing running spectrum analysis with chronux
% cmn 06/2011


%Event_Name_EEG='Lfpx'
Event_Name='pLFP'
max_events = 10^6;
max_t = 10^9;
LFPall=0;
for ch = 1:nChan
       
        [lfp ts] = readWave(TTX,ch, Event_Name,max_events,max_t);
        if timeSeriesOn
            lfpT{ch} = ts;
            lfpData{ch} = lfp;
        end
        if spectrumOn
        params.Fs = 1/median(diff(ts));
        params.tapers = [3 5];
        params.fpass = [0 150];
        
        [S t f] = mtspecgramc(lfp',[3 1],params);

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