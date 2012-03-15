%function readWaveforms
% Reads in waveform data and stores it into cluster file
% This is used to retroactively save out waveforms for data that has
% already been clustered, to take advantage of waveform comparison added to
% select_units by Erik Flister 03/12

close all
clear all
Tank_Name=[];

[fname pname] = uigetfile('','cluster file','*.mat');
load(fullfile(pname,fname));
pname
fname

timelimit=0;
maxT = 1000;

for tet=use_tets;
    clear X
    %%%% read in all waveforms for this tetrode
    for tet_ch=1:4
        N=0;
        Nblock=0;
        for block= 1:nblocks;
            
            ch = (tet-1)*4+tet_ch;
            data = getTDTdata(Tank_Name,char(Block_Name(block)),ch,flags);
            
            T = data.MUspikeT{ch};
            s = data.snips{ch};
            if timelimit
                T = T(T<maxT);
                s = s(:,1:length(T));
            end
            
            Nblock(block) = length(T);
            
            event_times_all{ch}(N+(1:Nblock(block))) = T +10^5 *(block-1);
            X(N+(1:Nblock(block)),1:size(s,1),tet_ch)=s';
            size(X);
            N=N+Nblock(block);
            N
            NS(tet_ch)= N;
        end
        
    end
    clear W
    noise = squeeze(std(X(:,1,:),[],1));
    
    %%% waveform data is now stored in X
    %%% first dimension is waveform #, second dimension is timepoint
    %%% third dimension is channel # (1-4)
    
    %%% sometimes TDT loses spikes.
    %%% if one channel is dropped for a snippet, then all the following
    %%% snippets are offset.
    %%% so we need to find any events that only have 3channels worth of
    %%% data, and remove them. we do this by looking at event times
    %%% note: shouldn't be necessary for a more reliable acquisition system!
    if std(NS)~=0
        %%% |  sum(std(event_times_all((tet-1)*4+1:4,:)))>0
        sprintf('warning missing spikes !!!')
        N=0
        
        Xold = X;
        event_times_old = event_times_all;
        
        spike_times= zeros(1,max(NS));
        X = zeros(size(X));
        for i = 1:NS(1)
            missing=0;
            match_t(1)=i;
            for test_ch = 2:4;
                m =find(event_times_old{(tet-1)*4+test_ch}==event_times_old{(tet-1)*4+1}(i));
                if ~isempty(m)
                    match_t(test_ch)=m;
                else
                    missing=1;
                end
            end
            if missing==0;
                N=N+1;
                for j=1:4
                    event_times_all{(tet-1)*4+j}(N) = event_times_old{(tet-1)*4+1}(i);
                    X(N,:,j) = Xold(match_t(j),:,j);
                end
                
            end
        end
        X=X(1:N,:,:);
        for i = 1:4
            event_times_all{(tet-1)*4+i} = event_times_all{(tet-1)*4+i}(1:N);
        end
    end
    
    
    
    %%% cut out beginning of waveform before threshold, no information
    %%% there (necessary?)
    X = X(:,6:30,:);
    
    clear Xshift Xold etimes_old newtimes;
    pack
    
    %%% begin alignment of spikewaveforms
    
    trigT = 6;  %%% timepoint to align on
    shiftrange=6;  %%% maximum shift allowed in alignement
    
    sz = 25-shiftrange;
    Xshift = zeros(N,sz,4);
    
    %%% remove snippets that have an initial voltage that is beyond the
    %%% usual noise range, as they are generally artifacts
    
    threshold_voltage = 1;  %%% use thresholding?
    v0thresh = 40*10^-6;    %%% threshold for initial point (a real spike will not have such a large deviation before beginning of spike
    vmax_thresh = 800*10^-6;  %%% threshold for largest excursion - only artifacts will have amplitude greater than this
    
    v0 = squeeze(max(abs(X(:,1:3,:)),[],3));
    v0 = squeeze(max(v0,[],2));
    vmax = squeeze(max(abs(X(:,:,:)),[],3));
    vmax = squeeze(max(vmax,[],2));
    
    if threshold_voltage
        used = find(v0< v0thresh & vmax <vmax_thresh);
    else
        used = find(v0>0);
    end
    figure
    hist(v0);
    
    
    %%% perform alignement based on minimum point, as averaged across channels
    for chan = 1:4
        [Y I] = min(squeeze(mean(X(:,trigT:trigT+shiftrange,:),3)),[],2);
        for i = 1:N
            Xshift(i,:,chan) = X(i,I(i) : (I(i)+sz-1),chan);
        end
    end
    
    wave_all{tet}=Xshift;
    size(Xshift);
    size(idx_all{(tet-1)*4+1})
    sprintf('finished tet %d ',tet)
end %tet
%
 save(fullfile(pname,fname))