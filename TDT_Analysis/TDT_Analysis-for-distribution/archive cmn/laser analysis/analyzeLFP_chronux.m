function [lfpT lfpData spectT spectF spectData] = analyzeLFP_chronux(TTX,nChan,timeSeriesOn,spectrumOn);
% Matlab code for reading LFP data and 
% performing running spectrum analysis with chronux
% cmn 06/2011


Event_Name_EEG='Lfpx'
max_events = 10^6;
max_t = 10^9;
LFPall=0;
for ch = 1:nChan
        N = invoke(TTX, 'ReadEventsV',    max_events , Event_Name_EEG,ch, 0, ...
            0,max_t,'All');
        if (N==max_events)
            warning('max number of events aquired');
        end
        W = invoke(TTX, 'ParseEvV', 0, N);
        Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        W=W(:);     
        Wave_TS = Wave_TS(:);
        lfpData{ch} = W';
        lfpT{ch} = Wave_TS;
        
        dt=mean(diff(Wave_TS))/32;
        
        tstamp = zeros(length(W),1);
        for i = 1:length(Wave_TS);
            tstamp((i-1)*32+1:i*32) = Wave_TS(i)+(0:31)*dt;
        end
        tstamp = tstamp(1:length(W));
        
        params.Fs = 1/dt;
        params.tapers = [3 5];
        params.fpass = [0 90];
        
        [S t f] = mtspecgramc(W,[3 1],params);
        
        Snorm = zeros(size(S));
        for i = 1:size(S,1);
            Snorm(i,:)=(S(i,:)).*((1:size(S,2)).^1);
        end
        df = max(f)/length(f);
        
        spectData{ch} = S;
        spectT{ch} = t;
        spectF{ch} = f;
        
        %     figure
        %     imagesc(Snorm',[0 10^-7]);
        %     axis xy
        %     set(gca,'YTick',(10:10:80)/df);
        %     set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
        
        
        figure
        %imagesc(Snorm,[0 0.1]);
        imagesc(Snorm',[0 6*10^-8]);
        axis xy
        set(gca,'YTick',(10:10:80)/df);
        set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    
end %% tet