%function compile_analysis_pkmzeta data
clear all
close all

[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


for dataset = 1:2
    
    cells=0;
    if dataset ==1
        afiles = {'C:\data\matlab data\Jen data_analysis\12_3_12\record_1\analysis_1.mat',...
            'C:\data\matlab data\Jen data_analysis\12_3_12\record_2\analysis_2.mat',...
            'C:\data\matlab data\Jen data_analysis\12_5_12\record_1\analysis_1a.mat',...
            'C:\data\matlab data\Jen data_analysis\12_5_12\record_2\analysis_2a.mat'};
    elseif dataset ==2
        afiles = {'C:\Users\nlab\Desktop\developmental_periods\12_10_12\analysis_1.mat',...
           'C:\Users\nlab\Desktop\developmental_periods\12_10_12\record_2\analysis_2.mat'}; %%% tg
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
            [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        end
    end
    
   % sprintf('mean firing rate = %f',mean(peak))
    meanpeak(dataset) = mean(peak)
    meanOSI(dataset) = mean(OSI(peak>2))
    meanwidth(dataset) = mean(width(peak>2 & OSI>0.5))
    
    meanspont(dataset) = mean(driftspont(peak>2))
    meanwpref(dataset) = mean(driftwpref(peak>2 &driftwpref<=0.32))
    medianwpref(dataset) = median(driftwpref(peak>2 &driftwpref<=0.32))
    meanF1F0(dataset) = mean(driftF1F0(peak>2));
   %to generate variables to do stats on
   bothwpref{dataset}=driftwpref(peak>2 &driftwpref<=0.32);
   bothpeak{dataset}=peak;
   bothOSI{dataset}=OSI(peak>2);
   bothwidth{dataset}=width(peak>2 &OSI>0.5);
   bothspont{dataset}=driftspont(peak>2);
   bothF1F0{dataset}=driftF1F0(peak>2);

    % to generate error bars for SEM on graphs
    peak_err(dataset) = std(peak)/sqrt(length(peak));
    OSI_err(dataset) = std(OSI(peak>2))/sqrt(length(OSI(peak>2)));
    width_err(dataset) = std(width(peak>2 &OSI>0.5))/sqrt(length(width(peak>2 &OSI>0.5)));
    wpref_err(dataset) = std(driftwpref(peak>2 &driftwpref<=0.32))/sqrt(length(driftwpref(peak>2 &driftwpref<=0.32)));
    spont_err(dataset) = std(driftspont(peak>2))/sqrt(length(driftspont(peak>2)));
    F1F0_err{dataset} = std(driftF1F0(peak>2))/sqrt(length(driftF1F0(peak>2)));
end

figure
barweb(meanpeak,peak_err);
title('Mean Peak Amplitude')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%

figure
barweb(meanOSI,OSI_err);
title('Mean OSI')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwidth,width_err);
title('Mean Width of Orientation Selective Response')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwpref,wpref_err);
title('Mean Preferred Spatial Frequency')

figure
barweb(medianwpref,wpref_err);
title('Median Preferred Spatial Frequency')

figure
barweb(meanspont,spont_err);
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
