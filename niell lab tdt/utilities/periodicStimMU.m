  pname = uigetdir('C:\data\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    nchan = input('number of channels : ');
    flags = struct('visStim',1,'MUspike',1,'snips',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    
    frame_duration = 1/60;
     
    
    close all
    
    clear n_spikes ntrials
    cyc_all=0;
    for ch = [53 61]
        sp = data.snips{ch};
        amp =min(sp);
        sptimes = data.MUspikeT{ch};
        thresh = prctile(amp,50);
        sprintf('thresh = %f',thresh*10^6)
        sptimes = sptimes(amp<thresh);
       % sptimes = sptimes(amp<thresh & amp>-150*10^-6);
        
    for f = 1:max(data.frameEpocs(1,:))
                [Spike_Timing index numtrials] = getTrialsSU(data.frameEpocs,sptimes, f, frame_duration);
            n_spikes(f) = length(Spike_Timing);
            ntrials(f) = numtrials;
    end        
             n_spikes = n_spikes./(frame_duration*ntrials);      
             ch
%              figure
%              plot(n_spikes)
             for f = 1:600;
                 cyc_avg(f) = mean(n_spikes(f+240:600:end));
             end
%              figure
%              plot(cyc_avg)
             clear cyc
             for t=1:50
                 cyc(t) = median(cyc_avg(12*(t-1)+1 :12*t));
             end
             figure
             plot(cyc)
             title(sprintf('%d',ch));
             
             cyc_all = cyc_all + cyc;
    end
    
   figure
   plot(cyc_all)