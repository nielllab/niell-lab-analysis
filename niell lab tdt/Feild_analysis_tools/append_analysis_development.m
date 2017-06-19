%assign waveform ID 
clear all
close all
dbstop if error

apath = 'F:\Jennifer_Development\'; %path on niellV2 %apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0; wvform_type=0; 
sessionNum=0;

for dataset = 1:2  %%% control ("wt") NR5A1-cre/CHR2 animals vs.NR2B vs. NR2A deleted NR5A1-cre
    
    if dataset ==1
        
        
        afiles = { 'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
           'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
            'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
            'Good recordings\Adults\9_12_13\rec1\analysis_9_12_13_adult_rec1.mat',...
            'Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',... %%% no wn
            'Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
            'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
            'Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
            'Good recordings\Adults\11_14_13\rec1\analysis_11_14_13_adult_rec1.mat',...
            'Good recordings\Adults\11_14_13\rec2\analysis_11_14_13_adult_rec2.mat',...
             'Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
             'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'};
 
    elseif dataset==2
%         
        afiles = { 'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',... 
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...           
            'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...    
            'Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec1\analysis_11_20_13_rec1.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec2\analysis_11_20_13_EO0_rec2.mat',...
            'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
            'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat' };

%     elseif dataset==3
%         afiles = {'N2B\8_19_15_rec2_analysis_2.mat',...
%            'N2B\analysis_12_16_15.mat' }
%      
    end
    
    for i = 1:length(afiles)
        
        clear trough2peak n_units cells trough_depth peak_height inh midnarrow exc OSI DSI peak width psth ...
            dpeak driftA1 driftA2 driftB drift_theta_w drift_theta driftwpref driftF1F0cvDSI cvOSI
       
        afile = fullfile(apath,afiles{i});
        load([apath afiles{i}]);
       
        afiles{i}
        
        n_units = length(cells);
       
        %%% determine waveform ID
%         F1 = trough2peak;
%         F2 = -trough_depth./peak_height;
%         EI=[F1;F2]'
%        
%         opts = statset('display','final');
%         [idx,ctrs] = kmeans(EI,2,'distance','city','Replicates',5,'options',opts);
%         if sum(idx==1)<sum(idx==2)
%             inh = (idx==1);
%         else
%             inh = (idx==2);
%         end

       
        trough2peak;
        TH_ratio = -trough_depth./peak_height;
        inh = TH_ratio>0 & TH_ratio <1.6  & trough2peak<9;  %%% directly based on wvform; k-means includes another inh group?
        %midnarrow = TH_ratio>0 & TH_ratio<4 & trough2peak<10 &  trough2peak(:,5)>7.5;  %%% could analyze these specifically at some point
        exc= trough2peak>9 & TH_ratio>1.5 & TH_ratio<4;
      
        figure
        plot(trough2peak,TH_ratio,'ko');
        hold on
        plot(trough2peak(find(inh)),TH_ratio(find(inh)),'ro');hold on
        plot(trough2peak(find(exc)),TH_ratio(find(exc)),'go');

        inh=inh';
        exc=exc';
        %determine whether unit is light responsive
        if  exist('drift', 'var'); 
        
        %dpeak=field2array(drift,'maxFR'); %determine whether unit is responisve to drifting gratings
        driftA1= field2array(drift,'A1');
        driftA2=field2array(drift,'A2');
        driftB= field2array(drift,'B');
        drift_theta_w=field2array(drift,'thetawidth');
        drift_theta=field2array(drift,'theta');
        driftwpref = field2array(drift,'wpref');
        layer=field2array(drift,'layer');
        driftF1F0= field2array(drift,'F1')./field2array(drift,'F0');
        
%         cvDSI = field2array(drift,'cv_dsi'); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
%         cvOSI=field2array(drift,'cv_osi');
        else
            
        %dpeak=NaN;
        driftA1= NaN;
        driftA2=NaN;
        driftB= NaN;
        drift_theta_w=NaN;
        drift_theta=NaN;
        driftwpref = NaN;
        driftF1F0= NaN;
        layer=NaN;
%         cvDSI = NaN; 
%         cvOSI=NaN;
       
        end
          
        
        for j = 1:size(driftA1,1);
            for k=1:size(driftA1,2);
                driftA1(j,k);
                driftA2(j,k);
                driftB(j,k);
                drift_theta_w(j,k);
                [OSI(j,k) DSI(j,k) width(j,k) peak(j,k)] = calculate_tuning(driftA1(j,k),driftA2(j,k),driftB(j,k),drift_theta_w(j,k));
                
                
            end
        end
        
        resp = peak(:,1)>=1;
        OS=OSI(:,1);
        %cvOS=cvOSI(:,1);
        layer=layer(:,1);
        TW=width(:,1);
        prefOrient=drift_theta(:,1);
        DS=DSI(:,1);
        %cvDS=cvDSI(:,1);
        SF=driftwpref(:,1);
        
      
        %load in peri-light stim firing rates
%         psth_pinp=psth;
%          %histbins = -50ms to + 50 in steps of 1 (msec) if you need to plot
%          
%          plotrange = 48:90;
%          figure
%          plot(plotrange,psth_pinp(:,plotrange),'k');
%         
%         baseline = mean(psth_pinp(:,5:45),2);
%         baseStd = std(psth_pinp(:,5:45),[],2);
% 
%         ev = max(psth_pinp(:,52:55),[],2);
%         evoked = ev- baseline;
% 
%         zscore =evoked./baseStd;
%         zscore(zscore>20)=20;
% 
%         pinped = (zscore>5& evoked>15 & ~inh);
        
        save (afile,'resp','layer', 'OS','TW','prefOrient','DS','SF','inh','-append')
    
       
    end %%% loop over adult vs EO
end