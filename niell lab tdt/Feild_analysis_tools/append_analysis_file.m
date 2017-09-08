%assign waveform ID and Whether ChR2 tagged
clear all
close all
dbstop if error

%cd 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\'
apath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\'; %apath = 'D:\Jen_ephys_data\developmental_periods\';

N =0; cells=0;PINPed=0; wvform_type=0; 
sessionNum=0;

for dataset = 1  %%% control ("wt") NR5A1-cre/CHR2 animals vs.NR2B vs. NR2A deleted NR5A1-cre
    
    if dataset ==1
        
        
        afiles = { '1_13_15_analysis_2.mat',...
        '1_18_16_analysis_2.mat',...   
        '2_25_15_analysis_2.mat',...
           '3_11_15_analysis_2.mat',...
          '4_9_15_analysis_2.mat',...
         '4_10_15_analysis_2.mat',...
         '4_13_15_analysis_2.mat',...
         '4_30_15_analysis_2.mat',...
         '5_4_15_analysis_2.mat',...
         '5_11_15_analysis_2.mat',...
         '6_25_15_analysis_2',...
         '6_29_15_analysis_2A',...
         '8_17_15_analysis_2.mat'};
 %'NR5A_Pinping\8_17_15\analysis_2.mat',... no pinped cell met criterion
 %'8_10_15_analysis_2.mat',...
 %'NR5A_Pinping\2_19_15\analysis.mat',...
 %'NR5A_Pinping\6_17_15\analysis_2',...
 %'NR5A_Pinping\4_18_15\analysis_2.mat',...
 %'NR5A_Pinping\3_13_15\analysis_2.mat',...
 %'NR5A_Pinping\3_24_15\full\analysis_2.mat',...
 %'NR5A_Pinping\3_7_15\analysis_2.mat',...
    elseif dataset==2
%         
        afiles = { 'N2A\3_3_15_analysis_2.mat',...
            'N2A\3_4_15_analysis_2.mat',...
            'N2A\4_24_15_analysis_2A.mat',...
            'N2A\5_12_15_analysis_2.mat',...
            'N2A\6_28_15_analysis_2A',...
            'N2A\6_30_15_analysis_2B',...
            'N2A\7_29_15_analysis_2',...
            'N2A\8_18_15_analysis_2.mat',...
            'N2A\8_18_15_rec2_analysis_2.mat',...
            'N2A\9_12_15_analysis.mat',...
            'N2A\9_16_15_analysis_2.mat',...
            'N2A\11_2_15_analysis_2.mat',...
            'N2A\11_3_15_analysis_2.mat',...
            'N2A\analysis_1_20_16.mat'};
        %'N2A\8_6_15_analysis_2.mat',...
%         %'N2A\8_7_15_analysis_2A.mat',...no pinped units
    elseif dataset==3
        afiles = { 'N2B\2_2_15_analysis_2.mat',...
            'N2B\2_24_15_analysis_2.mat',...
            'N2B\3_25_15_analysis_3_25_15.mat',...
            'N2B\3_26_15_analysis_2.mat',...
            'N2B\4_23_15_analysis_2.mat',...
            'N2B\6_18_15_analysis_2A.mat',...
            'N2B\6_22_15_analysis_2A.mat',...
            'N2B\6_23_15_analysis_2.mat',...
            'N2B\6_26_15_analysis_2A.mat',...
            'N2B\7_1_15_analysis_2.mat',...
            'N2B\8_11_15_analysis_2A.mat',...
            'N2B\8_11_15_rec2_analysis_2.mat',...
            'N2B\8_13_15_analysis_2A.mat',...
            'N2B\8_13_15_rec2_analysis_2A.mat',...
            'N2B\8_14_15_analysis.mat',...
            'N2B\8_19_15_analysis_2.mat',...
            'N2B\8_19_15_rec2_analysis_2.mat',...
           'N2B\analysis_12_16_15.mat' };
%      
    end
    
    for i = 1:length(afiles)
        
        afile = fullfile(apath,afiles{i});
        load([apath afiles{i}]);
       
        afiles{i}
        
    
        %%% determine waveform ID
        clear F1 F2 EI inh idx ctrs
        F1 = trough2peak;
        F2 = -trough_depth./peak_height;
        EI=[F1;F2]'
       
        opts = statset('display','final');
        [idx,ctrs] = kmeans(EI,2,'distance','city','Replicates',5,'options',opts);
        if sum(idx==1)<sum(idx==2)
            inh = (idx==1);
        else
            inh = (idx==2);
        end

       
%         trough2peak;
%         TH_ratio = -trough_depth./peak_height;
%         inh = TH_ratio>0 & TH_ratio <1.6  & trough2peak<9;  %%% directly based on wvform; k-means includes another inh group?
%         midnarrow = TH_ratio>0 & TH_ratio<4 & trough2peak<10 &  trough2peak(:,5)>7.5;  %%% could analyze these specifically at some point
%         exc= trough2peak>9 & TH_ratio>1.5 & TH_ratio<4;
      
%         figure
%         plot(trough2peak,TH_ratio,'ko');
%         hold on
%         plot(trough2peak(find(inh)),TH_ratio(find(inh)),'ro');hold on
%         plot(trough2peak(find(exc)),TH_ratio(find(exc)),'go');
        
        figure 
           plot(F1,F2,'ko'); hold on
           plot(F1(find(inh)),F2(find(inh)),'ro')
           

        

        inh_revised=inh';
        %exc_revised=exc';
        %determine whether unit is light responsive
        
        if  exist('drift', 'var'); 
        
        dpeak=cell2mat(field2array(drift,'maxFR')); %determine whether unit is responisve to drifting gratings
        driftA1= cell2mat(field2array(drift,'A1'));
        driftA2=cell2mat(field2array(drift,'A2'));
        driftB= cell2mat(field2array(drift,'B'));
        drift_theta_w=cell2mat(field2array(drift,'thetawidth'));
        drift_theta=cell2mat(field2array(drift,'theta'));
        driftwpref = cell2mat(field2array(drift,'wpref'));
        driftF1=cell2mat(field2array(drift,'F1'));
        driftF0=cell2mat(field2array(drift,'F0'));
        DriftF1F0=driftF1./driftF0;
        cvDSI = cell2mat(field2array(drift,'cv_dsi')); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
        cvOSI=cell2mat(field2array(drift,'cv_osi'));
        else
            
        dpeak=NaN;
        driftA1= NaN;
        driftA2=NaN;
        driftB= NaN;
        drift_theta_w=NaN;
        drift_theta=NaN;
        driftwpref = NaN;
        driftF1F0= NaN;
        cvDSI = NaN; 
        cvOSI=NaN;
       
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
        
        resp = peak(:,1)>=2;
        OS=OSI(:,1);
        cvOS=cvOSI(:,1);
        TW=width(:,1);
        prefOrient=drift_theta(:,1);
        DS=DSI(:,1);
        cvDS=cvDSI(:,1);
        SF=driftwpref(:,1);
        
       % load in peri-light stim firing rates
        psth_pinp=psth;
         %histbins = -50ms to + 50 in steps of 1 (msec) if you need to plot
         
         plotrange = 48:90;
         figure
         plot(plotrange,psth_pinp(:,plotrange),'k');
        
        baseline = mean(psth_pinp(:,5:45),2);
        baseStd = std(psth_pinp(:,5:45),[],2);

        ev = max(psth_pinp(:,52:54),[],2);
        evoked = ev- baseline;

        zscore =evoked./baseStd;
        zscore(zscore>20)=20;

        pinped_revised = (zscore>8& evoked>25 & ~inh);
        
       save (afile,'pinped_revised','resp', 'OS','cvOS','TW','prefOrient','DS','cvDS','SF','inh_revised','-append')
    
       end
  % end %%% loop over adult vs EO
end