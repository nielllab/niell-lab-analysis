%function compile_analysis_pkmzeta data
for dataset = 1:2
    
    cells=0;
    if dataset ==1
        afiles = {'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\12_22_11\1st_penetration_recluster\analysis.mat',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\12_22_11\2ndpen\analysis_2nd_penetration',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_04_12\bars16D_drift_orient_sf\analysis_bars_drift_orient_sf.mat',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_04_12\3rd_record\analysis_3rd.mat',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_09_12\2ndpenetration\analysis.mat',...
        'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_09_12\3rd_penetration_polytrode\bars_drift_3rd.mat',...
       'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_17_12\1st_record\select_.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_08_12\1st Penetration\analysis_1.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_08_12\2nd\analysis.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_08_12\3rd\abalysis_3.mat',...
     'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_28_12\3rd_recording\analysis_3rdRecording.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_28_12\4th_recording\analysis_4th.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\03_01_12\analysis_1.mat',...
     'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\03_01_12\1a\1A_analysis.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\03_01_12\2ndrecord\select_analysis.mat'}
    elseif dataset ==2
        afiles = {'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_16_12\1stpenetration\select_1.mat',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_16_12\2nd pentration\select_2.mat',...
        'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_16_12\4th penetration\select_4_analysis.mat',...
         'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_31_12\1st_pentration\bars_drift_WN_1\analysis_bars_WN_1.mat',...
        'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_31_12\2nd_penetration\drift_bars_WN_2_analysis.mat',...
       'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\01_31_12\3rd_penetration\drift_bars_Wn_3analysis.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_27_12\analysis_1.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_28_12\1st_recording\analysis_1.mat',...
      'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_28_12\2nd recording_new site\analysis_2ndrecording.mat',...
       'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_29_12\1st_record\analysis_1strecord.mat',...
       'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_29_12\2c_recording_3rd_recording\analysis_2c.mat',...
       'C:\Users\nlab\Desktop\Jen_Hoy_PKMzeta\02_29_12\2nd_record\analysis_2ndrecord.mat'}; %%% tg
    end
    
    N =0;
    for i = 1:length(afiles)
        
        i
        load(afiles{i});
        n_units = length(L_ratio);
        cellrange = N+1:N+n_units;
        N=N+n_units;
        
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
        
        %     for c = 1:n_units;
        %         ch = alldata(c,1);
        %         cl = alldata(c,2);
        %         min_t = min(mean_wvform (:,ch : ch+3,cl),[],1);
        %         [trough_allc trig_chan] = min(min_t);
        %         trig_chan = ch+trig_chan-1;
        %         alldata(c,7) = mean_wvform(size(mean_wvform,1),trig_chan,cl)/peak_height(c);
        %
        %         t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == cl)));
        %         dt = diff(t1);
        %         dt =dt(dt<.02);
        %         n=hist(dt,.001:0.002:.02);
        %         [y alldata(c,8)] = max(n);
        %         n=hist(dt,.0005:0.001:.02);
        %         alldata(c,9) = max(n(3:8))./mean(n(15:20))
        %     end;
        
        
        %         A1(cellrange,:)=bars_A1;
        %         A2(cellrange,:)=bars_A2;
        %         w(cellrange,:)=bars_w;
        %         theta(cellrange,:)=bars_theta;
        %         bspont(cellrange,:)=barspont;
        %
        %         size(rf_width)
        %
        %         rfw(cellrange,:) = rf_width*30;
        
        
        size(driftorient_A1)
        driftA1(cellrange,:) =driftorient_A1;
        driftA2(cellrange,:) = driftorient_A2;
        driftB(cellrange,:) = driftorient_B;
        drift_theta_w(cellrange,:)=driftorient_thetawidth;
        driftspont(cellrange,:) = drift_spont;
        
        driftwpref(cellrange,:) = wpref_dog;
        driftwbw(cellrange,:) = wbw_dog;
        
        driftF1F0(cellrange,:) = driftorient_F1./driftorient_F0;
        
        %size(wvform)
        wvform(cellrange,:) = wv';
        
    end
    
    for i = 1:N
        for j=1:1
            driftA1(i,j);
            driftA2(i,j);
            driftB(i,j);
            drift_theta_w(i,j);
            [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j))
        end
    end
    
    sprintf('mean firing rate = %f',mean(peak))
    meanpeak(dataset) = mean(peak)
    meanOSI(dataset) = mean(OSI(peak>2))
    meanwidth(dataset) = mean(width(peak>2 & OSI>0.5))
    meanspont(dataset) = mean(driftspont(peak>2))
end





