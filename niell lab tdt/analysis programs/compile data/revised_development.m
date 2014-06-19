%function compile_analysis_pkmzeta data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


for dataset = 1:2
    
    cells=0;
    if dataset ==1
        afiles = {'D:\Jen_ephys_data\FragileX\1_22_13\Record1_drift_1b and wn_1\analysis_1_22_13_1A.mat',...
            'D:\Jen_ephys_data\FragileX\1_22_13\Record_2_drift_2 and wn_2\analysis_2_1_22_13.mat',...
           'D:\Jen_ephys_data\FragileX\01_24_13\A\Record_2_ist good recording\analysis_2.mat'};
    elseif dataset ==2
        afiles = {'D:\Jen_ephys_data\developmental_periods\10_12_12_P14_EO1\1st_recording\analysis_1.mat',...
            'D:\Jen_ephys_data\developmental_periods\10_13_12_P15_EO2\1st_record\analysis_1.mat',...
            'D:\Jen_ephys_data\developmental_periods\10_13_12_P15_EO2\2ndrecord\analysis_10_13_12_EO2_cluster_record2.mat',...
            'D:\Jen_ephys_data\developmental_periods\10_15_12_EO4\record_1\amalysis_1_10_15_12_EO2_cluster_1.mat'}; %%% tg
    end
    
    N =0;
    for i = 1:length(afiles)
        
        i;
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
        
        
        size(drift)
        driftA1(cellrange,:) =field2array(drift,'A1');
      
        driftA2(cellrange,:) =field2array(drift,'A2') ;
        driftB(cellrange,:) = field2array(drift,'B');
        drift_theta_w(cellrange,:)=field2array(drift,'thetawidth');
        driftspont(cellrange,:) = field2array(drift,'spont');
        
        driftwpref(cellrange,:) = field2array(drift,'wpref');
        driftwbw(cellrange,:) = field2array(drift,'bw') ;
        
        driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
        
        %size(wvform)
        wvform(cellrange,:) = wv';
        
    end
    
    for i = 1:size(driftA1,1)
        for j=1:size(driftA1,2)
            driftA1(i,j);
            driftA2(i,j);
            driftB(i,j);
            drift_theta_w(i,j);
            [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        end
    end
    
   % sprintf('mean firing rate = %f',mean(peak))
   
   %meanpeak(dataset,1) be stationary and meanpeak(dataset,2) 
  
   responsive = peak(:,1)>2 & peak(:,2)>2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
   responsive_tuned = peak(:,1)>2 & peak(:,2)>2 & OSI(:,1)>0.5 & OSI (:,2)>0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
   responsive_tuned_spatial = peak(:,1)>2 & peak(:,2)>2 & driftwpref(:,1) <=0.32 & driftwpref(:,2) <=0.32
   
   for moving =1:2 % moving variable describe state and 1:2 describes that there are two states to be considered for this analysis
   meanpeak(dataset,moving) = mean(peak(:,moving));
   
   meanOSI(dataset,moving) = mean(OSI(responsive,moving))

   meanwidth(dataset,moving) = mean(width(responsive_tuned,moving))
   
   meanspont(dataset,moving) = mean(driftspont(responsive,moving))
   
   meanwpref(dataset,moving) = mean(driftwpref(responsive_tuned_spatial,moving)) 
   
   medianwpref(dataset,moving) = median(driftwpref(responsive_tuned_spatial,moving))
   
   meanF1F0(dataset,moving) = mean(driftF1F0(responsive,moving));
   
   %to generate variables to do stats on
   bothwpref{dataset,moving}=driftwpref(responsive,moving);
   bothpeak{dataset,moving}=peak(responsive,moving);
   bothOSI{dataset,moving}=OSI(responsive,moving);
   bothwidth{dataset,moving}=width(responsive_tuned,moving);
   bothspont{dataset,moving}=driftspont(responsive,moving);
   bothF1F0{dataset}=driftF1F0(responsive,moving);

    % to generate error bars for SEM on graphs
    peak_err(dataset,moving) = std(peak(:,moving))/sqrt(length(peak(:,moving)));
   
    OSI_err(dataset,moving) = std(OSI(responsive,moving))/sqrt(length(OSI(responsive,moving)));
    width_err(dataset,moving) = std(width(responsive_tuned,moving))/sqrt(length(width(responsive_tuned,moving)));
    wpref_err(dataset,moving) = std(driftwpref(responsive_tuned_spatial,moving))/sqrt(length(driftwpref(responsive_tuned_spatial,moving)));
    spont_err(dataset,moving) = std(driftspont(responsive,moving))/sqrt(length(driftspont(responsive,moving)));
    F1F0_err{dataset,moving} = std(driftF1F0(responsive,moving))/sqrt(length(driftF1F0(responsive,moving)));
  end
  
end

figure
barweb(meanpeak',peak_err');
title('Mean Peak Amplitude')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%

figure
barweb(meanOSI',OSI_err');
title('Mean OSI')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwidth',width_err');
title('Mean Width of Orientation Selective Response')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwpref',wpref_err');
title('Mean Preferred Spatial Frequency')

figure
barweb(medianwpref',wpref_err');
title('Median Preferred Spatial Frequency')

figure
barweb(meanspont',spont_err');
title('Mean Spontaneous Activity')

clear h
figure
histbins= 0:2:50;
for rep=1:2
    h(rep,:) = hist(bothspont{rep},histbins)/length(bothspont{rep});
end
bar(histbins,h')
legend('wt','tg')
title('histogram of spontaneous activity')
ylabel('firing rate') 

%figure
%histbins= -30:2:30;
%for rep=1:2
 %   h(rep,:) = hist(bothpeak{rep})/length(bothpeak{rep});
%end
%bar(histbins,h')
%legend('wt','tg')

clear h
figure
histbins= 0:0.1:1;
for rep=1:2
    h(rep,:) = hist(bothOSI{rep},histbins)/length(bothOSI{rep});
end
bar(histbins,h')
legend('wt','tg')
title('OSI')

%figure
%histbins= 0:0.2:1;
%for rep=1:2
 %   h(rep,:) = hist(bothwidth{rep},histbins)/length(bothwidth{rep});
%end
%bar(histbins,h')
%legend('wt','tg')


figure
histbins= [0.01 0.02 0.04 0.08 0.16 0.32];
for rep=1:2
    wh(rep,:) = hist(bothwpref{rep},histbins)/length(bothwpref{rep});
end
bar(1:6,wh')
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',histbins)
legend('wt','tg')
title('spatial frequency preference')

figure
histbins= linspace(0,2,10);
for rep=1:2
    linh(rep,:) = hist(bothF1F0{rep},histbins)/length(bothF1F0{rep});
end
bar(histbins,linh')
legend('wt','tg')
title('F1F0')

%%% convert ps to pdf and delete old ps
ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);


% figure
% barweb(meanF1F0,F1F0_err);
% title('Mean Periodicity')



ranksum(bothpeak{1},bothpeak{2})
ranksum(bothOSI{1},bothOSI{2})
ranksum(bothwpref{1},bothwpref{2})
ranksum(bothwidth{1},bothwidth{2})
ranksum(bothspont{1},bothspont{2})
ranksum(bothF1F0{1},bothF1F0{2})
