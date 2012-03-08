clear alldata
alldata(:,1:2) = cells;
alldata(:,3) = L_ratio;

%%% waveform
alldata(:,4) = trough_width;
alldata(:,5) = trough2peak;
alldata(:,6) = -trough_depth./peak_height;

for c = 1:size(alldata,1);
    ch = alldata(c,1);
    cl = alldata(c,2);
    min_t = min(mean_wvform (:,ch : ch+3,cl),[],1);   
      [trough_allc trig_chan] = min(min_t);
       trig_chan = ch+trig_chan-1;
       alldata(c,7) = mean_wvform(size(mean_wvform,1),trig_chan,cl)/peak_height(c);
    
     t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == cl)));
    dt = diff(t1);
    dt =dt(dt<.02);
    n=hist(dt,.001:0.002:.02);
    [y alldata(c,8)] = max(n);
     n=hist(dt,.0005:0.001:.02);
     alldata(c,9) = max(n(3:8))./mean(n(15:20))
       
       
end;

xlswrite('wvform',alldata);
alldata

clear alldata
%%% bars 16d
%%% responsiveness
prefamp = (bars_B' + bars_A1');
%alldata(:,1) = prefamp./(prefamp + bars_spont);
%%% orientation selectivity
alldata(:,1) = (prefamp-bars_null')./(prefamp + bars_null');
%%% preferred orientation
alldata(:,2) = bars_theta'*180/pi;
%%% orientation widht
alldata(:,3) = abs(bars_w)*180/pi;
%%% direction selectivity
alldata(:,4) = (bars_A1'-bars_A2')./(bars_A1'+bars_A2' + 2*bars_B');
%%% peak response
alldata(:,5) = prefamp;
%%% spontaneous
alldata(:,6) = bars_spont;
alldata(:,7) = rf_width*30;

xlswrite('bars',alldata);
alldata


%%% drift gratings
clear alldata
prefamp = (driftorient_A1' +driftorient_B');
%%% responsiveness
alldata(:,1) = prefamp ./(drift_spont + prefamp)
%%% orientation selectivity
%driftorient_null(driftorient_null<0) =0;
alldata(:,2) = (prefamp - driftorient_null') ./ (prefamp + driftorient_null');
%%% preferred orientation
alldata(:,3) = mod(driftorient_theta+pi,2*pi)*180/pi;
%%% orientation width (can be negative since gaussian is symmetric
alldata(:,4) = abs(driftorient_thetawidth)*180/pi;
%%% direction selectivity
alldata(:,5) = (driftorient_A1 - driftorient_A2)./(driftorient_A1 );
%%% spatial frequnecy
alldata(:,6) = driftorient_wpref;
%%% spatial frequency width
alldata(:,7) = driftorient_bw  %%% sigma in factors of two
%%% F1F0
alldata(:,8) = driftorient_F1 ./ driftorient_F0;
%%% peak response
alldata(:,9) = prefamp;
%%% spontaneous
alldata(:,10) = drift_spont;
alldata(:,11:12)=cells;
xlswrite('drift',alldata(:,2:12));
alldata

clear alldata
alldata(:,1:2) = cells;
alldata(:,3) = bars_OSI';
%alldata(:,4) = bars_dsi';
alldata(:,5) = driftorient_osi';
alldata(:,6) = driftorient_dsi';
xlswrite('circvar',alldata);

clear alldata;
%%%% movies
prefamp = (sta_A1' + sta_baseline');
%%% preferred orientation
alldata(:,1) = mod(sta_thetapref,pi)*180/pi;
%%% preferred SF
%%% fft is taken over 48 pix, minimum frequency = 1/48 pix
%%%60 pixels = 70 degrees; 48pix = 48*60/60 = 56 degrees
%alldata(:,2) = sta_wpref/(48*70/60); 
%%% note now it's taken over full 60x60 (with zero padding outside
alldata(:,2) = sta_wpref/(70); 
%%% total spikes
alldata(:,3) = sta_N;
%%% orientation selectivity
alldata(:,4) = (prefamp - sta_null')./(prefamp + sta_null');
%%% responsiveness to contrast varying
alldata(:,5) = sta_responsiveness;
alldata(:,6) = mod(sta_contrastphase,2*pi);
alldata(:,7) = duration_wn;
alldata(:,8) = halfcontrast_wn;
alldata(:,9) = halfslope_wn;
alldata(:,10) = adaptation_index;
xlswrite('sta',alldata);
alldata

clear alldata

%xlswrite('contrast',response_wn');



% 
%%% shortbar
clear alldata
alldata(:,1) =sb_x0mean;
alldata(:,2)=sb_y0mean;
alldata(:,3)=sb_wx0mean;
alldata(:,4) = sb_wy0mean;
alldata(:,5) = sb_ampmean;
alldata(:,6) = sb_spont;
xlswrite('sb',alldata);
alldata

%%% counterphase
clear alldata
alldata(:,1:4) = cptf_tftuning;
alldata(:,5) = cptf_F0;
alldata(:,6) = cptf_F1;
alldata(:,7) = cptf_F2;
alldata(:,8) = cptf_spont;
xlswrite('cphase',alldata);
alldata

%%% flash
clear alldata
alldata(:,6) =flash_used;
alldata(:,1) = flash_amp;
alldata(:,2)=flash_on;
alldata(:,3) = flash_off;
alldata(:,4) = flash_phasedependence;
alldata(:,5) = on_off_phasedependence;
xlswrite('flash',alldata);
alldata





clear alldata
alldata(:,1:2) = cells;
alldata(:,3) = drift_peakR;
alldata(:,4) = drift_peakvar;
alldata(:,5)=drift_peakerr;
alldata(:,6) = drift_spontvar;
alldata(:,7)= drift_sponterr;
alldata
xlswrite('drift_variance',alldata);

