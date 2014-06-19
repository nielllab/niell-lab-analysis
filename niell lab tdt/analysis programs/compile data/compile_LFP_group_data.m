%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;
for dataset = 1:2  %%% adult vs eye opening
    
    if dataset ==1
    afiles = {'Good recordings\EO7_EO9\8_6_13_EO7\analysis_EO7_rec1_8_6_13.mat',...
              'Good recordings\EO7_EO9\8_6_13_EO7\rec2\analysis_8_6_13_rec2_EO7.mat',...
              'Good recordings\EO7_EO9\9_2_13_EO7\rec2\analysis_9_2_13_EO9_rec2.mat',...
              'Good recordings\EO7_EO9\05_10_13_EO7\Record1_upper\analysis_1.mat',...
              'Good recordings\EO7_EO9\05_10_13_EO7\Record_1_deeper\analysis_EO8_deeper.mat',...
              'Good recordings\EO7_EO9\9_6_13_EO11\analysis_9_6_13_EO11_rec1.mat',...
              'Good recordings\EO7_EO9\9_6_13_EO11\rec2\analysis_9_6_13_EO11_rec2.mat',...
              'Good recordings\EO7_EO9\11_27_13_EO7\rec1\analysis_11_27_13_EO7_rec1.mat',...
              'Good recordings\EO7_EO9\11_27_13_EO7\rec2\analysis_11_27_13_EO7_rec2.mat'};
    elseif dataset ==2
        afiles = { 'Good recordings\EO3_EO4\5_6_13_EO3\mouseC\analysis_mouseC_EO3_5_6_13.mat',...
             'Good recordings\EO3_EO4\7_5_13P16_EO4\mouseB\analysis_07_5_13_cluster_7_5_13_mouseB_945uM_analysis.mat',...
             'Good recordings\EO3_EO4\8_2_13_EO3\analysis_EO3_8_2_13_rec1.mat',...
             'Good recordings\EO3_EO4\9_10_13_EO3\rec1\analysis_9_10_13_EO3_rec1.mat',...
             'Good recordings\EO3_EO4\9_10_13_EO3\rec2\analysis_9_10_13_EO3_rec2.mat',...
             'Good recordings\EO3_EO4\11_25_13_EO4\rec1\analysis_11_25_13_EO4_rec1.mat',...
             'Good recordings\EO3_EO4\11_25_13_EO4\rec2\analysis_11_25_13_EO4_rec2.mat'}; %% no wn
         
    end
    %   if dataset ==1
%        afiles = {'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
%            'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
%             'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
%             'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
%             'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
%             'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
%             'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
%             'Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',... %%% no wn
%             'Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
%             'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
%             'Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
%             'Good recordings\Adults\11_14_13\rec1\analysis_11_14_13_adult_rec1.mat',...
%             'Good recordings\Adults\11_14_13\rec2\analysis_11_14_13_adult_rec2.mat',...
%              'Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
%              'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'}; %% no wn
%      elseif dataset ==2
%          afiles = { 'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
%             'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',... 
%             'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...           
%             'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...    
%             'Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat',...
%             'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
%             'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat' }; %%% tg
%        end
    
%        'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
%             'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',... 
%             'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...           
%             'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...    
%             'Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat',...
%             'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
%             'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat',...  

%            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
%            'Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
%            'Good recordings\EO1_EO2\7_17_13_EO1\Rec2\Analysis_7_17_13_rec2_EO1.mat',...  
%            'Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
%            'Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
%            'Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
%            'Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
%            'Good recordings\EO1_EO2\11_20_13_EO0\rec1\analysis_11_20_13_rec1.mat',...
%            'Good recordings\EO1_EO2\11_20_13_EO0\rec2\analysis_11_20_13_EO0_rec2.mat',...
%            'Good recordings\EO1_EO2\11_27_13_mouseB_EO2\analysis_11_27_13_EO2_mouseB.mat'
       
    for i = 1:length(afiles)
        
        
        clear wn wn_movement
        clear LFP_movement
        
        load([apath afiles{i}]);
        n_units = length(L_ratio);
        cellrange = N+1:N+n_units;
        N=N+n_units;
        
        number(i) = n_units;
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
        alldata( cellrange,7:25)= wv';
        
        age(cellrange)=3-dataset;
        
        driftlayer =  field2array(drift,'layer');
        lyr(cellrange,:) = driftlayer(:,1);
        

        %size(wvform)
        wvform(cellrange,:) = wv';
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        
        
        ev=[]; sp=[];lfp = []; stopcrf=[];mvcrf=[];lfp_bar = [];
        if exist('wn','var')
            for j = 1:length(wn);
             
               if ~isempty(wn(j).N)
                    ev(j) = mean(wn(j).crf(9:12))-mean(wn(j).crf([1 2 19 20]));
                    sp(j) =mean(wn(j).crf([1 2 19 20]));
                    if exist('wn_movement','var')
                        for mv = 1:2
                        lfp(j,mv,:) = interp1(wn_movement(j).freqs,wn_movement(j).mv_lfp(mv,:),1:120);
                        end
                   stopcrf(j,:)=wn_movement(j).stopCRF;
                   mvcrf(j,:) = wn_movement(j).moveCRF;
                    else
                        lfp=NaN; stopcrf=NaN; mvcrf=NaN;
                        sprintf('no wn movement!!')
                    end
                
                else
                    
                    ev(j)=NaN;
                    sp(j)=NaN;
                    lfp(j,:,:)=NaN; mvcrf(j,:)=NaN;stopcrf(j,:)=NaN;
                    
                end
            end
            moveLFP(cellrange,:,:)=lfp;  
            wn_evoked(cellrange)=ev;
            wn_spont(cellrange)=sp;
            wn_mv(cellrange,:)=mvcrf;
            wn_stop(cellrange,:)=stopcrf;
        else
            display('no wn!!')
            afiles{i}

         wn_mv(cellrange,:)=NaN;
         wn_stop(cellrange,:)=NaN;

           moveLFP(cellrange,:,:)=NaN;
           wn_evoked(cellrange)=NaN;
           wn_spont(cellrange)=NaN;
         
         
        end
          
   if exist ('LFP_movement', 'var')
    
      for k = 1:length(LFP_movement)   
           if ~isempty(LFP_movement(k).mv_lfp) 
           for mv = 1:2
                lfp_bar(k,mv,:) = interp1(LFP_movement(k).freqs,LFP_movement(k).mv_lfp(mv,:),1:120);
           end 
           else 
                lfp_bar(k,:,:) = NaN;
                sprintf('no LFP_movement');
           end     
      end
  
   moveLFP_bar(cellrange,:,:)=lfp_bar;  
   
  else
       
   moveLFP_bar(cellrange,:,:)=NaN;
  
   end
  
    end %%loop over cells
         
end %%% loop over adult vs EO

        
    
% LFPmax = max(nanmax(moveLFP,[],3),[],2);
% figure
% hist(LFPmax)
%moveLFP(:,:,59:61)=0.5*moveLFP(:,:,59:61);

LFPmax_Bar = max(nanmax(moveLFP_bar,[],3),[],2);
figure
hist(LFPmax_Bar)
%moveLFP_bar(:,:,60:61)=0.5*moveLFP_bar(:,:,60:61);
moveLFP_bar_1=moveLFP_bar(:,:,1:80);
moveLFP_bar_1(:,:,60)=1.1*moveLFP_bar_1(:,:,58);
moveLFP_bar_1(:,:,59)=1.1*moveLFP_bar_1(:,:,57);
moveLFP_bar_1(:,:,61)=1.1*moveLFP_bar_1(:,:,62);
LFPmax_Bar = max(nanmax(moveLFP_bar_1,[],3),[],2);
figure
hist(LFPmax_Bar)

%layerAgePlotMv([wn_stop(:,1) wn_mv(:,1)],age,lyr,inh,1,'spont');

%%% peak firing rate in gratings by age and movement)
%layerAgeActivity(peak(:,1),peak(:,2),age,lyr,inh,1,{'stationary','moving'},'drift evoked');


% figure
% hist(peak(:,2),-40:40);hold on
% xlabel('drift peak moving');
% 
% figure
% hist(peak(:,1),-40:40);hold on
% xlabel('drift peak stat');

%%% spont firing rate for gratings by age and movement
% layerAgeActivity(driftspont(:,1),driftspont(:,2),age,lyr,inh,1,{'stationary','moving'},'drift spont');
% 
% 
% [ m e n] =layerAgePlotMv([wn_stop(:,10)-wn_stop(:,1),wn_mv(:,10)-wn_mv(:,1)],age,lyr,inh,1,'evoked');
% 
% [ m e n] =layerAgePlotMv([peak(:,1),peak(:,2)],age,lyr,inh,peak(:,1)>1 & peak(:,2)>1,'drift evoked');
% 
% [ m e n] =layerAgePlotMv([driftspont(:,1),driftspont(:,2)],age,lyr,inh,1,'drift spont');

%plot the LFP power spectrum for the whole population during WN stim
% figure
% imagesc(squeeze(moveLFP(:,1,:)))
% 
% figure
% imagesc(squeeze(moveLFP(:,2,:)))
% 
% figure
% imagesc(squeeze(moveLFP(LFPmax<1*10^4,2,:)))
% 
% figure
% imagesc(squeeze(moveLFP(LFPmax<1*10^4,1,:)))

%%bar stim LFPs

% figure
% imagesc(squeeze(moveLFP_bar(:,1,:)))

lfp=LFPmax_Bar<10.0*10^4 & LFPmax_Bar>100;

%%%all units, normalized to largest signal




figure
moveLFP_bar_all_run=squeeze(moveLFP_bar_1(lfp & age'==2& lyr==5,2,:));
moveLFP_bar_all_run=sort(moveLFP_bar_all_run);
imagesc(moveLFP_bar_all_run,[0 prctile(moveLFP_bar_all_run(:),92)]);
title ('all units running')

figure
moveLFP_bar_all_stat=squeeze(moveLFP_bar_1(lfp & age'==2 & lyr==5,1,:));
moveLFP_bar_all_stat =sort(moveLFP_bar_all_stat);
imagesc(moveLFP_bar_all_stat,[0 prctile(moveLFP_bar_all_stat(:),92)]);
title ('all units stationary')


figure
moveLFP_bar_all_run=squeeze(moveLFP_bar_1(lfp & age'==1 & lyr<=3,2,:));
[~ ,ind_moveLFP]=sortrows(moveLFP_bar_all_run,[57 32 7]);
moveLFP_bar_all_run=moveLFP_bar_all_run(ind_moveLFP,:)
imagesc(moveLFP_bar_all_run,[0 prctile(moveLFP_bar_all_run(:),92)]);
title ('all units running')

figure
moveLFP_bar_all_stat=squeeze(moveLFP_bar_1(lfp & age'==1 & lyr<=3,1,:));
moveLFP_bar_all_stat =moveLFP_bar_all_stat(ind_moveLFP,:);
imagesc(moveLFP_bar_all_stat,[0 prctile(moveLFP_bar_all_stat(:),92)]);
title ('all units stationary')


%%%make compilation LFP power figs based on layer,age,stat info
%moveLFP_bar=squeeze(moveLFP_bar(LFPmax_Bar<1.0*10^4 & age'==1 & lyr<=6,2,:));


%%%running data EO1

%plot the LFP power running vs. stationary for EO1 and adult populations 


figure
for ly = 2:6
    plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==1 & lyr<=ly & lfp,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==1& lyr<=ly& lfp,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median EO1')


%%%EO1




for ly =2:6
    figure
    N=sum(~isnan(moveLFP_bar_1(age'==1 & lyr==ly & lfp)));
    mean_lfp= squeeze(nanmedian(moveLFP_bar_1(age'==1 & lyr==ly & lfp,1,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==1 & lyr==ly & lfp,1,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'b');
   % errorbar(mean_lfp,SEM_lfp,'b')
   % plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:),1)),'Linewidth',1+ly/3);
hold on
    mean_lfp= squeeze(nanmedian(moveLFP_bar_1(age'==1 & lyr==ly & lfp,2,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==1 & lyr==ly & lfp,2,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'g');
%plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median EO3 L6')

figure
for ly =3:6
    figure
    N=sum(~isnan(moveLFP_bar_1(age'==1 & lyr==ly & lfp)));
    mean_lfp= squeeze(nanmean(moveLFP_bar_1(age'==1 & lyr==ly & lfp,1,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==1 & lyr==ly & lfp,1,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'r');
   % errorbar(mean_lfp,SEM_lfp,'b')
   % plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:),1)),'Linewidth',1+ly/3);
hold on
    mean_lfp= squeeze(nanmean(moveLFP_bar_1(age'==1 & lyr==ly & lfp,2,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==1 & lyr==ly & lfp,2,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'k');
%plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean EO3_L6')
%%%LFP bars pop avg
%plot the LFP power running vs. stationary for EO1 and adult populations 


figure

for ly =3:6
    figure
    N=sum(~isnan(moveLFP_bar_1(age'==2 & lyr==ly & lfp)));
    mean_lfp= squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'k');
   % errorbar(mean_lfp,SEM_lfp,'b')
   % plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:),1)),'Linewidth',1+ly/3);
hold on
    mean_lfp= squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'r');
%plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median adult L3')

figure
for ly = 3
    figure
    N=sum(~isnan(moveLFP_bar_1(age'==2 & lyr==ly & lfp)));
    mean_lfp= squeeze(nanmean(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==1 & lyr==ly & lfp,1,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'r');
   % errorbar(mean_lfp,SEM_lfp,'b')
   % plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,1,:),1)),'Linewidth',1+ly/3);
hold on
    mean_lfp= squeeze(nanmean(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:)));
    stddv=squeeze(nanstd(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:)));
    SEM_lfp=stddv/(sqrt(N));
    shadedErrorBar(1:80,mean_lfp,SEM_lfp,'k');
%plot(1:80,squeeze(nanmedian(moveLFP_bar_1(age'==2 & lyr==ly & lfp,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean EO1_L3')

