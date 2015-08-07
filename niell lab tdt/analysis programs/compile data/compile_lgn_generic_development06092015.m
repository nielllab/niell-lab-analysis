%function compile_lgn_fulldataset
close all
clear all
n_obs=0;
day = [1 1 1 2 3 3 ];

PostnatalAge = [17 17 16 18 18 22 16 17 17 17 22 22 60 60 16 16 18 18 19 19 60 14 14 14 19 19 20 20 18 23 25 25 17 17 18 18 19 19 60 60 60 60 60]
PostnatalAge(PostnatalAge==60)=30;
figure
hist(PostnatalAge,min(PostnatalAge):max(PostnatalAge))

figure
hist(PostnatalAge(PostnatalAge<60),14:25)

agelist=sort(unique(PostnatalAge));
agelist=agelist(agelist~=25);


analysisPath ='E:\Wayne Matlab data\Wayne Matlab data_fromC\analysis files\development lgn\'
afile = { '05022013_analysis_recording_1_P17.mat'...
    '05022013_analysis_recording_2_P17.mat'...
    '05072013_analysis_recording_1_P16.mat'...
    '05092013_analysis_recording_1_P18.mat'...
    '05092013_analysis_recording_2_P18.mat'...
    '05132013_analysis_recording_1_P22.mat'...
    '05232013_analysis_recording_1_P16.mat'...
    '05242013_analysis_recording_1_P17.mat'...
    '05242013_analysis_recording_2_P17.mat'...
    '05242013_analysis_recording_3_P17.mat'...
    '05282013_analysis_recording_1_P22.mat'...
    '05282013_analysis_recording_2_P22.mat'...
    '06052013_analysis_recording_1_adult.mat'...
    '06052013_analysis_recording_2_adult.mat'...
    '07052013_analysis_recording_1_P16.mat'...
    '07052013_analysis_recording_2_P16.mat'...
    '07072013_analysis_recording_1_P18'...
    '07072013_analysis_recording_2_P18.mat'...
    '07082013_analysis_recording_1_P19'...
    '07082013_analysis_recording_2_P19.mat'...
    '07152015_analysis_recording_1_adult.mat'...
    '07182013_analysis_recording_2_P14.mat'...
    '0718a2013_analysis_recording_1_P14.mat'...
    '0718a2013_analysis_recording_2_P14.mat'...
    '08122013_analysis_recording_1_P19.mat'...
    '08122013_analysis_recording_2_P19.mat'...
    '08132013_analysis_recording_1_P20.mat'...
    '08132013_analysis_recording_2_P20.mat'...
    '09052013_analysis_recording_2_P18.mat'...
    '09102013_analysis_recording_1_P23.mat'...
    '09122013_analysis_recording_1_P25.mat'...
    '09122013_analysis_recording_2_P25.mat'...
    '09242013_analysis_recording_1_P17.mat'...
    '09242013_analysis_recording_2_P17.mat'...
    '09252013_analysis_recording_1_P18.mat'...
    '09252013_analysis_recording_2_P18.mat'...
    '09262013_analysis_recording_1_P19.mat'...
    '09262013_analysis_recording_2_P19.mat'...
    '03122013_analysis_recording_2_adult.mat'...
    '0307a2013_analysis_recording_1_adult.mat'...
    '12072012_analysis_recording_1_adult.mat'...
    '02242015_analysis_rec1_adult.mat'...
    '02252015_analysis_rec1_adult.mat'...
    };


histoPath = 'C:\data\matlab data\Wayne Matlab data\histology\development histology\'
histfile = { '','';...
    '','';...
    '05072013_P16_sec2_sl2_recording_1_LGNsection3_electrode1to32.tif','05072013_P16_sec2_sl2_recording_1_LGNsection3_electrode33to64.tif';  ...
    '05092013_P18_sec1slide2_recording_1_LGN section3_electrode1to32.tif','05092013_P18_sec1slide2_recording_1_LGN section3_electrode33to64.tif';...
    '05092013_P18_sec2slide2_recording_2_LGN section4_electrode1to32.tif','05092013_P18_sec2slide2_recording_2_LGN section4_electrode33to64.tif';...
    '','';... 
    '','';...
    '','';...
    '','';...
    '','';...
    '','';...
    '','';...
    '06052013_adult_sec3slide2_recording_1_LGN section_electrode1to32.tif','06052013_adult_sec3slide2_recording_1_LGN section_electrode33to64.tif';...
    '','';...
    '','';...
    '','';...
    '07072013_P18_sec5slide1_recording_1_LGN section_electrode1to32.tif','07072013_P18_sec5slide1_recording_1_LGN section_electrode33to64.tif';...
    '07072013_P18_sec5slide1_recording_2_LGN section_electrode1to32.tif','07072013_P18_sec5slide1_recording_2_LGN section_electrode33to64.tif';...
    '07082013_P19_sec5slide2_recording_1_LGN section_electrode1to32.tif','07082013_P19_sec5slide2_recording_1_LGN section_electrode33to64.tif';...
    '07082013_P19_sec4slide2_recording_2_LGN section_electrode1to32.tif','07082013_P19_sec4slide2_recording_2_LGN section_electrode33to64.tif';...
    '','';...
    '','';...
    '','';...
    '','';...
    '08122013_P17_sec6slide1_recording_1_LGN section_electrode1to32.tif','08122013_P17_sec6slide1_recording_1_LGN section_electrode33to64.tif';...
    '08122013_P17_sec6slide1_recording_2_LGN section_electrode1to32.tif','08122013_P17_sec6slide1_recording_2_LGN section_electrode33to64.tif';...
    '08132013_P18_sec2slide2_recording_1_LGN section_electrode1to32.tif','08132013_P18_sec2slide2_recording_1_LGN section_electrode33to64.tif';...
    '08132013_P18_sec1slide2_recording_2_LGN section_electrode1to32.tif','08132013_P18_sec1slide2_recording_2_LGN section_electrode33to64.tif';...
    '09052013_P18_sec7slide1_recording_2_LGN section5_electrode1to32.tif','09052013_P18_sec7slide1_recording_2_LGN section5_electrode33to64.tif';...
    '09102013_P23_sec5slide1_recording_1_LGN section4_electrode1to32.tif','09102013_P23_sec5slide1_recording_1_LGN section4_electrode33to64.tif';...
    '09122013_P25_sec1slide2_recording_1_LGN section4_electrode1to32.tif','09122013_P25_sec1slide2_recording_1_LGN section4_electrode33to64.tif';...
    '09122013_P25_sec2slide2_recording_2_LGN section5_electrode1to32.tif','09122013_P25_sec2slide2_recording_2_LGN section5_electrode33to64.tif';...
    '09242013_P17_sec2slide1_recording_1_LGN section2_electrode1to32.tif','09242013_P17_sec2slide1_recording_1_LGN section2_electrode33to64.tif';...
    '09242013_P17_sec4slide1_recording_2_LGN section4_electrode1to32.tif','09242013_P17_sec4slide1_recording_2_LGN section4_electrode33to64.tif';...
    '09252013_P18_sec4slide1_recording_2_LGN section2_electrode1to32.tif','09252013_P18_sec4slide1_recording_2_LGN section2_electrode33to64.tif';...
    '09252013_P18_sec5slide1_recording_1_LGN section3_electrode1to32.tif','09252013_P18_sec5slide1_recording_1_LGN section3_electrode33to64.tif';...
    '09262013_P19_sec1slide2_recording_1_LGN section4_electrode1to32.tif','09262013_P19_sec1slide2_recording_1_LGN section4_electrode33to64.tif';...
    '09262013_P19_sec1slide2_recording_2_LGN section4_electrode1to32.tif','09262013_P19_sec1slide2_recording_2_LGN section4_electrode33to64.tif';...
    '','';...
    '','';...
    '','';...
    '','';...
    '','';...
    };


for i = 1:length(afile)
     afile{i} = [analysisPath afile{i}];
    for s = 1:2  
    if ~isempty(histfile{i,s})
        histfile{i,s} = [histoPath histfile{i,s}];
    end
    end
end



%%% WWT
%%% note - this assumes both shanks are in roughly same A/P position

lgnPos  = [0 0 3 3 4 0 0 0 0 0 0 0 3 0 0 0 2 3 1 3 0 0 0 0 2 1 4 3 1 1 2 2 1 4 2 4 3 2 2 0 0];

%%% lgnPos = ceil(5*rand(length(afile),1));

length(afile)
length(histfile)
length(lgnPos)

[anatomy sections]=LGN_histo_twoshank(histfile, lgnPos);


numsites = 32*ones(1,length(histfile));


% manual_type = xlsread('C:\data\lgn rf project_new0824\lgn_analysis\lgn types 091412.xlsx','A1:A294');
% manual_outside = xlsread('C:\data\lgn rf project_new0824\lgn_analysis\lgn types 091412.xlsx','B1:B294');
% manual_outside(isnan(manual_outside))=0;
n=0;

for i = 1:length(afile);
    i
    clear displayOffset drift mv fl wn wn_movement
    load(afile{i});
    cell_range = n+(1:length(drift));
    n=n+length(drift);
    size(wn)
    if size(wn,2)==2
        wn_all(cell_range,:)=wn
    else
        wn_all(cell_range,1)=wn
    end
    drift_all(cell_range)=drift;
    if exist('mv','var')
        mv_all(cell_range)=mv;
    else
        mv_all(cell_range)=NaN;
    end
    fl_all(cell_range)=fl;
    if exist('wn_movement','var')
        if ~isfield(wn_movement,'spikes')
            wn_movement(1).freqs = NaN; wn_movement(1).mv_lfp = NaN; wn_movement(1).mvlfp_tdata=NaN;
            wn_movement(1).speed = NaN; wn_movement(1).spikes = NaN; wn_movement(1).lfpV = NaN; wn_movement(1).lfpT=NaN;
            wn_movement(1).frameNum=NaN; wn_movement(1).frameT=NaN;
        end
        movement_all(cell_range)=wn_movement;
    end
    
    age(cell_range) = PostnatalAge(i);
    site(cell_range)=i;
    cell_id(cell_range,:) = cells;
    if size(wv,1)>19
        wv=wv(1:19,:);
    end
    wvform(cell_range,:) = wv';
    if ~exist('displayOffset')
        display(sprintf('need to add display for file %d',i))
       % input('press return to continue')
        displayOffset =45;
        displayHeight=4;
    end
    offsetX(cell_range)=displayOffset;
    offsetY(cell_range)=displayHeight*atand(1/25);
    if numsites(i)==16;
        peakchan= peakchan+16;
    end
    peaksite(cell_range)=peakchan;
    clear x y z s
    for j =1:length(peakchan)
        if isempty(anatomy(i).AP)
            x(j) = NaN;
            y(j) = NaN;
            z(j)=NaN;
            s(j)=NaN;
        else
            x(j)=anatomy(i).siteXY(1,peakchan(j));
            y(j)=anatomy(i).siteXY(2,peakchan(j));
            z(j)=anatomy(i).AP;
            s(j)=anatomy(i).section;
        end
    end
    
    histox(cell_range)=x;
    histoy(cell_range)=y;
    histoz(cell_range)=z;
    histSection(cell_range)=s;
    
end

%%% align points on trace to reference
alignedX= nan(size(histox));
alignedY= nan(size(histox));
for cell_n=1:length(histox)
    if ~isnan(histox(cell_n))
        
        paxxy = sections(histSection(cell_n)).coords;
        tracexy = anatomy(site(cell_n)).LGN;  %%% coords are reversed
        widthX = max(tracexy(:,2))-min(tracexy(:,2));
        normX(cell_n) = (histox(cell_n)-min(tracexy(:,2))) / widthX;
        ypts = find(abs((tracexy(:,2))-(histox(cell_n)))<10);
        if length(ypts)>1
            normY(cell_n) = (histoy(cell_n)-min(tracexy(ypts,1))) / (max(tracexy(ypts,1)) - min(tracexy(ypts,1)));
        else
            normY(cell_n)=0;
        end
        alignedX(cell_n) = min(paxxy(:,1)) + normX(cell_n)*(max(paxxy(:,1))-min(paxxy(:,1)));
        ypts = find(abs((paxxy(:,1))-(alignedX(cell_n)))<10);
        if length(ypts)>1
            alignedY(cell_n) = min(paxxy(ypts,2))+normY(cell_n)*(max(paxxy(ypts,2))-min(paxxy(ypts,2)));
        else
            alignedY(cell_n)=0;
        end
    end
end

originalX=histox;
originalY=histoy;
histox=alignedX;
histoy=alignedY;

%%% find inside points
inside= zeros(n,1);
for cell_n= 1:n
    if ~isnan(histox(cell_n))
        
        sec = sections(histSection(cell_n)).coords;
        y= sec(find(round(sec(:,1))==round(histox(cell_n))),2);
        if ~(isempty(y) || histoy(cell_n)>max(y) || histoy(cell_n)<min(y))
            inside(cell_n)=1;
        end
    end
end

labels = zeros(n,3);
labels(:,3)=1;
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);

labels= ones(n,3);
labels(~inside,2)=0;
labels(~inside,3)=0;
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);


if ~exist('manual_outside','var')
    manual_outside = 1-inside;
end


labels= ones(n,3);
labels(~manual_outside,2)=0;
labels(~manual_outside,3)=0;
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);

if ~exist('manual_type','var')
    manual_type = ones(size(histox));
end

% %%% plot manual type on normal and aligned position
% for i = [0 1]
%     %for i = 0
%     id = nan(size(manual_type));
%     id(~isnan(manual_type))=0;
%     id(manual_type==i)=1;
%
%     [labels cmap clim]= makeColors(id,nan,'qual','Set1');
%     f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
%     colormap(cmap); set(gca,'Clim',clim); title(sprintf('type %d',i));
% end

wn_degperpix = wn(end).degperpix;

clear p
n_sta=0;
%%% fit wn movie stas
wn_cr=nan(n,2);
for eye=1:2;
    for cell_n = 1:n
        cell_n
        wn_all(cell_n,eye).sta_fit=[nan nan nan nan nan nan ];
        wn_all(cell_n,eye).fitsta=nan;
        if ~isempty(wn_all(cell_n,eye).N)
            
            wn_cr(cell_n,eye) = (wn_all(cell_n,eye).crf(10) - wn_all(cell_n,eye).crf(1));
            eye
            for s =1:3
                sta = double(squeeze(wn_all(cell_n,eye).svd_xy(s,:,:)));
                
                [fit g]= fitLGNrf(sta);
                background = find(abs(g-fit(2))<(0.1*abs(fit(1))));
                if fit(1)>1 | fit(3)<1 | fit(4)<1 | fit(3)>size(sta,1) | fit(4)>size(sta,2)
                    break
                end
                z=abs(fit(1))/std(sta(background))
                
                if z>5
                    wn_all(cell_n,eye).sta_t = wn_all(cell_n,eye).svd_t(:,s);
                    if wn_all(cell_n,eye).sta_t(6)<0;
                        wn_all(cell_n,eye).sta_t = wn_all(cell_n,eye).sta_t *-1;
                        fit(1) = fit(1)*-1;
                        fit(2)=fit(2)*-1;
                    end
                    wn_all(cell_n,eye).sta_fit=fit;
                    wn_all(cell_n,eye).sta_final=sta;
                    wn_all(cell_n,eye).zscore=z;
                    wn_all(cell_n,eye).fitsta= g;
                    %                 figure
                    %                 subplot(2,2,1)
                    %                 imagesc(sta); axis equal; axis tight;
                    %
                    %                 subplot(2,2,2)
                    %                 imagesc(g);axis equal; axis tight;
                    %
                    %                 subplot(2,2,4);
                    %                 plot(wn_all(cell_n,eye).sta_t)
                    %                 xlim([0 30])
                    %
                    
                    n_sta=n_sta+1;
                    break
                end
            end
        end
    end
end



for cell_n=1:n
wn_spont(cell_n) = mean(wn_all(cell_n,1).crf([1 20]))/600;
end
figure
hist(wn_spont);
% figure
% plot(dr_spont,wn_spont,'o')

% figure
% plot(burst_fraction(:,1),sf_amp,'o')

close all


%%% set age bins
ageBins = [14 16; 18 20; 22 24; 25 30]';
[sp sperr] = sortbyage(wn_spont',age,ageBins,1);
figure
errorbar(mean(ageBins,1),sp,sperr);
ylabel('wn spont rate (sp/sec)')

[sp sperr] = sortbyage(wn_spont',age,agelist,1);
figure
errorbar(agelist,sp,sperr);
ylabel('wn spont rate (sp/sec)')

[wnamp wnerr] = sortbyage(wn_cr(:,1)/600,age,ageBins,1);
figure
errorbar(mean(ageBins,1),wnamp,wnerr)
ylabel('evoked wn firing (sp/sec)')

[wnamp wnerr] = sortbyage(wn_cr(:,1)/600,age,agelist,1);
figure
errorbar(agelist,wnamp,wnerr)
ylabel('evoked wn firing (sp/sec)')


ind=0;
tmin=zeros(1,cell_n); tmax = tmin; wn_lat=tmin;

x0 = nan(n,1); y0 = nan(n,1); wx=nan(n,1),wy=nan(n,1);A = nan(n,1); wn_resp = nan(n,1); wn_phase=nan(n,1);
tmin= nan(n,1); tmax=nan(n,1); wn_lat = nan(n,1); wn_cr_norm=nan(n,1);wn_adapt=nan(n,1); dom_eye=ones(n,1)
mv_x0 = nan(n,1); mv_y0=nan(n,1); fl_x0=nan(n,1); fl_y0=nan(n,1);
dom_eye(wn_cr(:,2)>1000 & wn_cr(:,2)>wn_cr(:,1)& wn_cr(:,1)>-10^3)=2;
figure
scatter(wn_cr(:,1),wn_cr(:,2),6,dom_eye)
xlabel('contra eye response')
ylabel('ipsi eye response');
legend('contra unit','ipsi unit')

for cell_n = 1:length(histox)
    if ~isempty(movement_all(cell_n).moveCRF)
    move_crf(cell_n,:) = movement_all(cell_n).moveCRF;
    move_ev(cell_n) = (move_crf(cell_n,end) - move_crf(cell_n,1));
    stop_crf(cell_n,:) = movement_all(cell_n).stopCRF;
    stop_ev(cell_n) = (stop_crf(cell_n,end) - stop_crf(cell_n,1));
    burst_fraction(cell_n,:) = movement_all(cell_n).burst;
    spont_ratio(cell_n) = (move_crf(cell_n,1)-stop_crf(cell_n,1)) / (move_crf(cell_n,1)+stop_crf(cell_n,1));
    evoke_ratio(cell_n) = (move_ev(cell_n) - stop_ev(cell_n))/ (move_ev(cell_n) + stop_ev(cell_n));
    burst_ratio(cell_n) = (movement_all(cell_n).burst(2) - movement_all(cell_n).burst(1))/ (movement_all(cell_n).burst(2) + movement_all(cell_n).burst(1));
    else
     move_crf(cell_n,:)=NaN; move_gain(cell_n)=NaN; stop_crf(cell_n,:)=NaN; stop_gain(cell_n)=NaN;burst_fraction(cell_n,:)=NaN; burst_ratio(cell_n)=NaN;
    end
%     if length(movement_all(cell_n).spikes)>1
%         clear params
% params.Fs = 0.25*1/median(diff(movement_all(cell_n).lfpT(1:1000)));
%         %params.tapers = [50*movement_all(cell_n).lfpT(end) 5];
%         params.tapers = [5 10  1];
%         params.fpass = [0 120];
%         [C ph s12 s1 s2 t fre] = cohgramcpt(movement_all(cell_n).lfpV(1:4:end)',movement_all(cell_n).spikes',[10 10],params);
%        df = median(diff(fre)); dt = median(diff(t));
%        figure
%         subplot(2,2,1);
%         imagesc(C'); 
%         axis xy; set(gca,'Ytick',(0:20:90)/df); set(gca,'Yticklabel',num2str((0:20:90)'));
%         set(gca,'Xtick',(0:300:max(t))/dt);set(gca,'Xticklabel',num2str((0:5:max(t)/60)'));
%         xlabel('min'); ylabel('Hz');
%         
%         subplot(2,2,2);
%         imagesc(s1',[0 prctile(s1(:),90)]);
%          axis xy; set(gca,'Ytick',(0:20:90)/df); set(gca,'Yticklabel',num2str((0:20:90)'));
%         set(gca,'Xtick',(0:300:max(t))/dt);set(gca,'Xticklabel',num2str((0:5:max(t)/60)'));
%         xlabel('min'); ylabel('Hz');
%        
%         subplot(2,2,3)
%         imagesc(s2',[0 prctile(s2(:),95)]);
%          axis xy; set(gca,'Ytick',(0:20:90)/df); set(gca,'Yticklabel',num2str((0:20:90)'));
%         set(gca,'Xtick',(0:300:max(t))/dt);set(gca,'Xticklabel',num2str((0:5:max(t)/60)'));
%         xlabel('min'); ylabel('Hz');
%         
%         title(sprintf('cell %d',cell_n));
%         subplot(2,2,4)
%         plot(movement_all(cell_n).speed);
%         %         figure
% %         plot(fre, spikeLFPcoh(cell_n,:));
% %         figure
% %         plot(fre,s1)
% %             figure
% %         plot(fre,s2)
%     end   
end
age(age==60)=30;

[burst bursterror] = sortbyage(burst_fraction,age,agelist,~isnan(burst_fraction(:,1)') );
figure
errorbar(agelist(2:end),burst(2:end,1),bursterror(2:end,1),'r','Linewidth',2);  %%% double-check movement for p14
hold on
errorbar(agelist(2:end),burst(2:end,2),bursterror(2:end,2),'g','Linewidth',2);
legend({'stationary','moving'})

[burst bursterror] = sortbyage(burst_fraction,age,ageBins,~isnan(burst_fraction(:,1)') );
figure
errorbar(mean(ageBins,1),burst(:,1),bursterror(:,1),'r','Linewidth',2);  %%% double-check movement for p14
hold on
errorbar(mean(ageBins,1),burst(:,2),bursterror(:,2),'g','Linewidth',2);
legend({'stationary','moving'})


figure
hist(burst_fraction);
xlabel('burst fraction')

figure
hist(burst_ratio)
xlabel('burst_ratio');

figure
plot(burst_fraction(:,1),burst_fraction(:,2),'o');
xlabel('burst fraction stop')
ylabel('burst fraction move');
hold on
plot([0 0.5],[0 0.5]);

figure
plot(stop_crf(:,1),move_crf(:,1),'o');
xlabel('crf spont - stop')
ylabel('crf spont - move');
hold on
plot([0 10],[0 10])

figure
plot(stop_ev,move_ev,'o');
xlabel('crf evoked - stop');
ylabel('crf evoked - move');
hold on
plot([-10 10],[-10 10])



% gain used to be ev-spont/ ev + spont
% figure
% plot(stop_gain,move_gain,'o');
% xlabel('stop gain');
% ylabel('move gain');
% hold on
% plot([-1 1],[-1 1]);




%%% fill in white noise params
for cell_n=1:n
    for eye=1:2
        
        if ~isempty(wn_all(cell_n,eye).N)
            if eye==dom_eye(cell_n);
                %                 if eye==2
                %                     keyboard
                %                 end;
                dom_eye(cell_n)=eye;
                wn_resp(cell_n) = wn_all(cell_n,eye).responsiveness(1);
                wn_phase(cell_n) =wn_all(cell_n,eye).phase;
                A(cell_n) = wn_all(cell_n,eye).sta_fit(1);
                wx(cell_n)=abs(wn_all(cell_n,eye).sta_fit(6));
                wy(cell_n)=abs(wn_all(cell_n,eye).sta_fit(5));
                x0(cell_n)= (wn_all(cell_n,eye).sta_fit(4) - 64)*wn_degperpix + offsetX(cell_n);
                y0(cell_n)= (wn_all(cell_n,eye).sta_fit(3) - 38)*wn_degperpix ;  %%% top of monitor = -38, bottom = 38,  so moving upward means more negative
                y0(cell_n) = -y0(cell_n) + offsetY(cell_n);  %%% invert so y goes from bottom up
                if ~isempty(wn_all(cell_n,eye).sta_t)
                    [z wn_lat(cell_n)] = max(wn_all(cell_n,eye).sta_t);
                    tmax(cell_n)=max(wn_all(cell_n,eye).sta_t - wn_all(cell_n,eye).sta_t(1))*sign(A(cell_n));
                    tmin(cell_n)=min(wn_all(cell_n,eye).sta_t- wn_all(cell_n,eye).sta_t(1))*sign(A(cell_n));
                end
                %wn_cr(cell_n) = (wn_all(cell_n,eye).crf(10) - wn_all(cell_n,eye).crf(1));
                wn_cr_dom(cell_n) = wn_cr(cell_n,eye);
                wn_cr_norm(cell_n) = (wn_all(cell_n,eye).crf(10) - wn_all(cell_n,eye).crf(1)) ./ ( (wn_all(cell_n,eye).crf(10) + wn_all(cell_n,eye).crf(1)));
                up = mean(wn_all(cell_n,eye).crf(1:10)); dn = mean(wn_all(cell_n,eye).crf(11:20));
                if up~=0
                    wn_adapt(cell_n) = dn/up;
                end
                
                if ~isstruct(mv_all(cell_n))
                    mv_x0(cell_n)=NaN; mv_y(cell_n)=NaN;
                else
                    if ~isempty(mv_all(cell_n).sta_pos)
                        mv_x0(cell_n) = mv_all(cell_n).sta_pos(1)*mv_all(length(wn_all)).degperpix + offsetX(cell_n);
                        mv_y0(cell_n) = -mv_all(cell_n).sta_pos(2)*mv_all(length(wn_all)).degperpix +offsetY(cell_n);
                    end
                end
                
                if ~isstruct(fl_all(cell_n))
                    fl_x0(cell_n) = NaN; fl_y0(cell_n)=NaN;
                else
                    if ~isempty(fl_all(cell_n).sta_pos)
                        fl_x0(cell_n) = fl_all(cell_n).sta_pos(1)*fl_all(length(wn_all)).degperpix + offsetX(cell_n);
                        fl_y0(cell_n) = -fl_all(cell_n).sta_pos(2)*fl_all(length(wn_all)).degperpix +offsetY(cell_n);
                    end
                end
                
            end
        end
    end
end

% 
% figure
% plot (histox,x0,'o');
% xlabel('histology X');
% ylabel('RF X');
% 
% 
% figure
% plot (histoy,y0,'o');
% xlabel('histology Y');
% ylabel('RF Y');
% 
% labels=ones(n,3);
% goodfit = find(~isnan(x0));
% labels(goodfit,1)=0;
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% title('wn sta rfs');
% 
% [labels cmap clim]= makeColors(x0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('x0');
% 
% [labels]= makeColors(y0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('y0');


% [labels]= makeColors(mv_x0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('mv_x0');
%
% [labels]= makeColors(mv_y0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('mv_y0');
%
% [labels]= makeColors(fl_x0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('fl_x0');
%
% [labels]= makeColors(fl_y0,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% colormap(cmap); colorbar; set(gca,'Clim',clim); title('fl_y0');

%
% figure
% scatter(x0,y0,8, manual_type);
% xlabel('wn x0'); ylabel('wn y0')
%
% figure
% scatter(mv_x0,mv_y0,8, manual_type);
% xlabel('mv x0'); ylabel('mv y0')

% tratio=-tmin./tmax;
% tratio(tratio>1)=1;
% labels= makeColors(tratio,nan,'div','RdBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% title('min max ratio');

%%% flash spots
for cell_n = 1:n
    f= fl_all(cell_n);
    if ~isempty(f.N)
        flash_spont = mean(f.onset_hist(f.onset_bins<0.25))
        %     figure
        %      for r = 1:2
        %          hold on
        %          if r==1
        %          plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)))
        %          else
        %            plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'r')
        %          end
        %      end
        
        [junk fl_onoff(cell_n) ] = max(nanmean(f.hist_all,2));
        
        sz= squeeze(f.hist_all(fl_onoff(cell_n),:))-f.spont;
        fl_sztune(cell_n,:)=sz;
        [fl_amp(cell_n) fl_sz(cell_n)] = max(sz) ;
        fl_supp(cell_n) = sz(end)/max(sz);
        fl_spont(cell_n) = f.spont;
        c=corrcoef(squeeze(mean(f.onset_hist(1,:,:),2)), squeeze(mean(f.onset_hist(2,:,:),2)));
        fl_onoffcorr(cell_n) =c(2,1);
        use_sizes = find(sz>0.5*max(sz))
        for r=1:2
            
            h = squeeze(mean(f.onset_hist(r,use_sizes,:),2))-f.spont;
            fl_tseries(cell_n,r,:) = h;
            [onset(cell_n,r) onset_lat(cell_n,r)]=max(h(f.onset_bins>=0.3 & f.onset_bins<0.55));
            onset_sust(cell_n,r) = mean(h(f.onset_bins>=0.3 & f.onset_bins<0.55))/onset(cell_n,r);
            [offset(cell_n,r) offset_lat(cell_n,r)]=max(h(f.onset_bins>=0.55 & f.onset_bins<0.8));
            offset_sust(cell_n,r) = mean(h(f.onset_bins>=0.55 & f.onset_bins<0.8))/offset(cell_n,r);
        end
        
        
        if f.lag==1  %%% responds to offset
            fl_lat(cell_n)=offset_lat(cell_n,fl_onoff(cell_n));
            sustain(cell_n) =offset_sust(cell_n,fl_onoff(cell_n));
        else
            sustain(cell_n) =onset_sust(cell_n,fl_onoff(cell_n));
            fl_lat(cell_n)=onset_lat(cell_n,fl_onoff(cell_n));
        end
        fl_onset_amp(cell_n) = onset(cell_n,fl_onoff(cell_n));
        fl_lag(cell_n) = f.lag;
        fl_onoff(cell_n) = 2*(fl_onoff(cell_n)-1.5);
        if fl_lag(cell_n)==1  %%% if a cell responds with 1 frame delay, it is responding to offset
            fl_onoff(cell_n)=-1*fl_onoff(cell_n);
        end
    end
end

fl_thresh=1;

%%% size suppression
[supp supperr] = sortbyage(1-fl_supp',age,agelist,fl_amp>fl_thresh );
figure
errorbar(agelist,supp, supperr);
ylabel('size suppression');

[supp supperr] = sortbyage(1-fl_supp',age,ageBins,fl_amp>fl_thresh );
figure
errorbar(mean(ageBins,1),supp, supperr);
ylabel('size suppression');

%%% preferred spot size
[sz szerr] = sortbyage(fl_sz',age,agelist,fl_amp>fl_thresh);
figure
errorbar(agelist,sz, szerr);
ylabel('pref spot size')
 set(gca,'Ytick',1:6);
 set(gca,'Yticklabel',{'2','4','8','16','32','full'});

[sz szerr] = sortbyage(fl_sz',age,ageBins,fl_amp>fl_thresh);
figure
errorbar(mean(ageBins,1),sz, szerr);
ylabel('pref spot size')
 set(gca,'Ytick',1:6);
 set(gca,'Yticklabel',{'2','4','8','16','32','full'});


% [lat laterr] = sortbyage(fl_lat',age,agelist,fl_amp>fl_thresh);
% figure
% errorbar(agelist,lat, laterr);
% ylabel('flash spot latency (frames)')

%%% RF width
[resp resperr] = sortbyage(wx,age,agelist,1)
figure
errorbar(agelist,resp, resperr);
ylabel('wx from wn (sp/sec)')

[resp resperr] = sortbyage(wx,age,ageBins,1)
figure
hold on
errorbar(mean(ageBins,1),resp, resperr);
bar(mean(ageBins,1),resp);
ylabel('wx from wn (sp/sec)')

%%% spot reponse amplitude
[resp resperr] = sortbyage(fl_amp',age,agelist,1)
figure
errorbar(agelist,resp, resperr);
ylabel('spot response (sp/sec)')

[resp resperr] = sortbyage(fl_amp',age,ageBins,1)
figure
errorbar(mean(ageBins,1),resp, resperr);
ylabel('spot response (sp/sec)')

%%% fraction responsive to spots
[resp resperr] = sortbyage(fl_amp'>fl_thresh,age,agelist,1)
figure
errorbar(agelist,resp, resperr);
ylabel('% responsive flashing spots')

[resp resperr] = sortbyage(fl_amp'>fl_thresh,age,ageBins,1)
figure
errorbar(mean(ageBins,1),resp, resperr);
ylabel('% responsive flashing spots')



% fl_type= nan(size(fl_lag));
% fl_type(fl_lag==0)=1;
% fl_type(fl_lag==1)=2;
% labels= makeColors(fl_type,nan);
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% title('flash type')

%n_obs=n_obs+1; obs_name{n_obs} = 'fl size supp'; obs(:,n_obs) = fl_supp;


fl_sz(fl_amp<fl_thresh)=nan;
fl_lag(fl_amp<fl_thresh)=nan;
fl_type(fl_amp<fl_thresh)=nan;
sustain(fl_amp<fl_thresh)=nan;
fl_onoff(fl_amp<fl_thresh)=nan;
fl_lat(fl_amp<fl_thresh)=nan;
fl_onoffcorr(fl_amp<fl_thresh)=nan;

fl_x0(fl_amp<fl_thresh)=nan;
fl_x0(fl_x0<-20)=nan;
fl_y0(fl_amp<fl_thresh)=nan;

sustain(sustain<0)=nan;
sustain(sustain>1)=nan;

% labels= makeColors(fl_type,nan,'seq','GnBu');
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% title('fl_type')
% 
% sustain_norm = sustain;
% sustain_norm(sustain_norm<0)=0;
% sustain_norm(fl_onset_amp<2)=nan;
% labels= makeColors(sustain_norm);
% f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
% title('normalized sustain')

clear driftOSI driftDSI
mv_osi=nan(size(histox)); mv_prefO=nan(size(histox));
mv_dsi=nan(size(histox)); mv_prefD=nan(size(histox));
mv_spd = nan(size(histox));
mv_sz= nan(size(histox)); mv_spd_ratio=nan(size(histox));
mv_amp = nan(size(histox)); mv_sz_tune= nan(length(histox),3);
mv_sp_tune=nan(length(histox),5); mv_os_tune=nan(length(histox),8);
%%% drifting gratings & moving spots
for cell_n=1:length(x0)
    %for cell_n = [24]
    if isstruct(mv_all(cell_n))
        m= mv_all(cell_n);
        if ~isempty(m.N)
            [dr_amp dr_onoff(cell_n)] =max(mean(nanmean(nanmean(m.hist_all,4),3),2))
            onoff = fl_onoff(cell_n)*0.5 +1.5;
            onoff = dr_onoff(cell_n);
            
            h=squeeze(m.hist_all(dr_onoff(cell_n),:,:,:));
            %     if size(h,1)==5;
            %         h= h(2:4,:,:);
            %     end
            [x y z] = meshgrid(1:size(h,2),1:size(h,1),1:8);
            use = find(~isnan(h));
            missing = find(isnan(h));
            missingdata = griddatan([x(use) y(use) z(use)],h(use),[x(missing) y(missing) z(missing)]);
            filled_h = h;
            filled_h(missing) = missingdata;
            
            sz = squeeze(mean(nanmean(filled_h,3),2))-m.spont;
            use_sz = find(sz>0.5*max(sz));
            sp = squeeze(mean(nanmean(filled_h,3),1))-m.spont;
            use_sp = find(sp>0.5*max(sp));
            os_tune = squeeze(mean(nanmean(filled_h(use_sz,use_sp,:),2),1))-m.spont
            os_N= squeeze(sum(nansum(m.n_unique(onoff,use_sz,use_sp,:),3),2));
            
            mv_sz_tune(cell_n,:)=sz;
            mv_sp_tune(cell_n,:)=sp;
            mv_os_tune(cell_n,:)=os_tune;
            
            [mv_amp(cell_n) mv_sz(cell_n)] =max(sz);
            if length(sz)==5
                mv_sz(cell_n)=mv_sz(cell_n)-1;
            end
            [empt mv_spd(cell_n)] = max(sp);
            mv_spd_ratio(cell_n) = (sp(5)-sp(1))/(sp(5)+sp(1));
            [mv_osi(cell_n) mv_prefO(cell_n)]= calcOSI(os_tune,0);
            [mv_dsi(cell_n) mv_prefD(cell_n)]= calcOSI(os_tune,1);
            
            %         figure('Name',sprintf('cell %d',cell_n))
            %         subplot(2,2,1)
            %         plot(sz);
            %         axis([1 3 min(0,min(sz)) max(1,max(sz))])
            %         subplot(2,2,2)
            %         plot(sp);
            %         axis([1 5 min(0,min(sp)) max(1,max(sp))])
            %         subplot(2,2,3)
            %         plot(os_tune)
            %         axis([1 8 min(0,min(os_tune)) max(1,max(os_tune))])
            %         title(sprintf('OSI = %f DSI=%f',mv_osi(cell_n),mv_dsi(cell_n)))
            
        end
    end   
end

%%% drifting gratings
for cell_n=1:length(histox)
    d = drift_all(cell_n);
    dr_spont(cell_n)=d.spont(1);
    tfs = squeeze(nanmean(d.orient_tune(1:2,:,:),3));
    tf_ratio(cell_n,:) = (tfs(2,:)-tfs(1,:))./(tfs(2,:) + tfs(1,:));
    %figure('Name',sprintf('cell %d',cell_n))
    color = {'b','r'};
    [sf_amp(cell_n) peak_sf(cell_n)] = max(squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1)));
    
    drift_amp(cell_n) = max(max(max(d.orient_tune)));
    for r=1:2  %%% tempfreq
        lowsf(cell_n,r)=d.sf_tune(r,1,1);
        for f=1:2  %%% F0 F1
            [peak_amp peak_tf(cell_n)] = max(squeeze(nanmean(d.sf_tune(:,f,:),3)));
            drift_os = squeeze(d.orient_tune(r,f,:));
            [dr_osi(cell_n,r,f) dr_prefO(cell_n,r,f)]= calcOSI(drift_os,0);
            [dr_dsi(cell_n,r,f) dr_prefD(cell_n,r,f)]= calcOSI(drift_os,1);
            
            %                 subplot(2,2,f)
            %                 hold on
            %                 plot(drift_os,color{r});
            %                 title(sprintf('tfratio %f',tf_ratio(cell_n,f)));
            
            drift_sf= squeeze(d.sf_tune(r,f,:));
            %                 subplot(2,2,f+2)
            %                 hold on
            %                 plot(drift_sf,color{r});
            %                 title(sprintf('peak sf %d',peak_sf(cell_n)));
            %axis([1 8 min(0,min(drift_os)) max(drift_os)])
        end
    end
    pref_tf = round(sign(tf_ratio(cell_n,1))/2 + 1.5);
    if isnan(pref_tf)
        pref_tf = 1;
    end
    driftOSI(cell_n,:) = dr_osi(cell_n,pref_tf,:);
    driftDSI(cell_n,:) = dr_dsi(cell_n,pref_tf,:);
    driftOSItheta(cell_n,:) = dr_prefO(cell_n,pref_tf,:);
    driftDSItheta(cell_n,:) = dr_prefD(cell_n,pref_tf,:);
    driftF1F0(cell_n)=sum(d.orient_tune(pref_tf,2,:)) / sum(d.orient_tune(pref_tf,1,:)); %%% replace with max?
    %         driftF1F0(cell_n)=max(d.orient_tune(pref_tf,2,:)) / max(d.orient_tune(pref_tf,1,:)); %%% replace with max?
    lowpass(cell_n) = lowsf(cell_n,pref_tf);
end

%%% grating response amp
[resp resperr] = sortbyage(sf_amp',age,agelist,1);
figure
errorbar(agelist,resp,resperr);
ylabel('grating response (sp/sec)');

[resp resperr] = sortbyage(sf_amp',age,ageBins,1);
figure
errorbar(mean(ageBins,1),resp,resperr);
ylabel('grating response (sp/sec)');

%%% grrating responsive fraction
[resp resperr] = sortbyage(sf_amp'>1,age,agelist,1);
figure
errorbar(agelist,resp,resperr);
ylabel('% responsive gratings (>1sp/sec)');

[resp resperr] = sortbyage(sf_amp'>1,age,ageBins,1);
figure
errorbar(mean(ageBins,1),resp,resperr);
ylabel('% responsive gratings (>1sp/sec)');

%%% peak SF
sfs = [0 0.01 0.02 0.04 0.08 0.16 0.32 nan];
sf_inds = peak_sf;
sf_inds(peak_sf ==0) = 7;
sf_inds(isnan(peak_sf))=8;

[sf sferr] = sortbyage(sfs(sf_inds)',age,agelist,sf_amp>2)
figure
errorbar(agelist,sf,sferr);
ylabel('peak sf (cpd)');

[sf sferr] = sortbyage(sfs(sf_inds)',age,ageBins,sf_amp>2)
figure
errorbar(mean(ageBins,1),sf,sferr);
ylabel('peak sf (cpd)');

%%% fraction lowpass (greatest response to 0cpd)
[sf sferr] = sortbyage((peak_sf==1)',age,agelist,sf_amp>2)
figure
errorbar(agelist,sf,sferr);
ylabel('fraction lowpass');

[sf sferr] = sortbyage((peak_sf==1)',age,ageBins,sf_amp>2)
figure
errorbar(mean(ageBins,1),sf,sferr);
ylabel('fraction lowpass');

%%% gratings spont
[sp sperr] = sortbyage(dr_spont',age,agelist,1)
figure
errorbar(agelist,sp,sperr);
ylabel('spont rate gratings (sp/sec)');

[sp sperr] = sortbyage(dr_spont',age,ageBins,1)
figure
errorbar(mean(ageBins,1),sp,sperr);
ylabel('spont rate gratings(sp/sec)');

%%% fraction OS
[osi osierr] = sortbyage(driftOSI(:,1)>0.2,age,agelist,sf_amp>2);
figure
errorbar(agelist,osi,osierr);
ylabel('fraction orientation selective >0.2')

[osi osierr] = sortbyage(driftOSI(:,1)>0.2,age,ageBins,sf_amp>2);
figure
errorbar(mean(ageBins,1),osi,osierr);
ylabel('fraction orientation selective >0.2')

%%% fraction DS
[dsi dsierr] = sortbyage(driftDSI(:,1)>0.2,age,agelist,sf_amp>2);
figure
errorbar(agelist,dsi,dsierr);
ylabel('fraction direction selective >0.2')

[dsi dsierr] = sortbyage(driftDSI(:,1)>0.2,age,ageBins,sf_amp>2);
figure
errorbar(mean(ageBins,1),dsi,dsierr);
ylabel('fraction direction selective >0.2')

%%% fraction SBC
sbc= (wn_cr_dom<-1*10^3);

[sbc_frac sbc_err] = sortbyage(sbc',age,agelist,1);
figure
errorbar(agelist,sbc_frac,sbc_err);
ylabel('fraction sbc');

[sbc_frac sbc_err] = sortbyage(sbc',age,ageBins,1);
figure
errorbar(mean(ageBins,1),sbc_frac,sbc_err);
ylabel('fraction sbc');


keyboard

thresh=3;
mv_osi(mv_amp<thresh)=nan;
mv_dsi(mv_amp<thresh)=nan;
mv_spd(mv_amp<thresh)=nan;
mv_spd_ratio(mv_amp<thresh)=nan;
mv_x0(mv_amp<thresh)=nan;
mv_y0(mv_amp<thresh)=nan;

d_thresh=1;
driftF1F0(driftF1F0<0)=0;
tf_ratio(tf_ratio>1)=1;
tf_ratio(tf_ratio<-1)=-1;
tf_ratio(sf_amp<d_thresh)=nan;
peak_sf(sf_amp<d_thresh)=nan;
driftDSI(sf_amp<d_thresh,:)=nan;
driftOSI(sf_amp<d_thresh,:)=nan;
driftF1F0(sf_amp<d_thresh)=nan;




%%% lots of misc figs
for misc_figs = 1:0
    figure
    hist(wn_cr);
    
    figure
    hist(wn_cr_norm)
    
    figure
    polar(wn_phase,wn_resp,'o');
    figure
    plot(wn_resp.*cos(wn_phase),wn_cr_norm,'o');
    
    %%% wn amplitudes
    figure
    hist(A)
    sum(~isnan(A))/length(A)
    
    %%% wn sizes
    figure
    plot(wx*0.7,wy*0.7,'o');
    axis([0 10 0 10])
    axis square
    
    
    
    
    %%% wn on/off
    figure
    plot(tmax,tmin,'o')
    axis([-1 1 -1 1])
    axis square
    
    %%% wn latency
    figure
    hist(wn_lat(wn_lat>0))
    
    %%% wn sites
    figure
    color = {'bo','go','ro','co','mo','yo','ko'};
    for i=1:length(afile);
        hold on
        plot(x0(site==i),y0(site==i),color{mod(i,7)+1})
        
    end
    axis ij
    
    
    figure
    plot(tmin./tmax,wx,'o')
    
    figure
    plot(fl_lag,wx,'o')
    xlabel('flash_lag');
    ylabel('width deg');
    
    figure
    plot(x0,wx,'o')
    axis([0 128 0 15])
    xlabel('x - deg');
    ylabel('width x - deg');
    
    figure
    plot(x0,wy,'o')
    
    figure
    hist(fl_lat(fl_amp>4 & fl_lag==0));
    
    figure
    hist(sustain(fl_amp>4 & fl_lag==0));
    
    
    figure
    hist(tmin./tmax)
    
    figure
    plot(fl_lat(fl_amp>2 & fl_lag==0),sustain(fl_amp>2 & fl_lag==0),'o' );
    
    
    figure
    plot(tmax(fl_lag==0),tmin(fl_lag==0),'o')
    hold on
    plot(tmax(fl_lag==1),tmin(fl_lag==1),'ro')
    
    figure
    plot(fl_sz(fl_amp>4),wx(fl_amp>4),'o');
    hold on
    plot(fl_sz(fl_amp>4),wy(fl_amp>4),'go');
    figure
    plot(onset);
    figure
    plot(offset);
    
    figure
    plot(fl_onoff(fl_amp>4),tmax(fl_amp>4),'o')
    
    figure
    plot(onset(:,1),onset(:,2),'o')
    xlabel('off response onset')
    ylabel('on response onset')
    axis equal
    figure
    plot(onset(:,1),offset(:,2),'o');
    axis equal
    hold on
    plot(onset(:,2),offset(:,1),'ro')
    axis equal
    
    figure
    for i =1:2
        hold on
        if i==1
            plot(x0(fl_lag==0),y0(fl_lag==0),'o');
        elseif i==2
            plot(x0(fl_lag==1),y0(fl_lag==1),'go');
        end
    end
    title('position - lag =0 vs 1')
    axis ij
    
    fl_responder=find(onset>5)
    figure
    plot(onset_sust(fl_responder),onset_lat(fl_responder),'o')
    
    figure
    hist(onset_sust(onset>5));
    
    figure
    hist(onset_lat(onset>5),0:10)
    
    figure
    plot(fl_sz(fl_amp>1),fl_amp(fl_amp>1),'o')
    
    figure
    hist(fl_supp(fl_amp>2))
    
    figure
    plot(fl_lag(fl_amp>4),fl_supp(fl_amp>4),'o')
    
    
    fl_resp = find(fl_amp>4 | (fl_amp>2 & fl_spont<1))
    
end

for histo_analysis=1:0
    
    fl_y0 = fl_y0 - (nanmean(fl_y0)-nanmean(y0));
    mv_y0 = mv_y0 - (nanmean(mv_y0)-nanmean(y0));
    fl_x0 = fl_x0 - (nanmean(fl_x0)-nanmean(x0));
    mv_x0 = mv_x0 - (nanmean(mv_x0)-nanmean(x0));
    
    meanx0 = nanmean([x0 fl_x0 mv_x0],2);
    meany0 = nanmean([y0 fl_y0 mv_y0],2);
    
    offsets = zeros(length(x0),max(day));
    for i = 1:max(day);
        offsets(:,i)=(day(site)==i)';
    end
    
    %xyz = [ones(size(x0)) alignedX' alignedY' histSection'*200];
    xyz = [offsets alignedX' alignedY' histSection'*200];
    
    xc = regress(meanx0,xyz);
    yc = regress(meany0,xyz);
    xfit = xyz*xc;
    yfit = xyz*yc;
    
    colordef white
    figure
    scatter(meanx0,xfit,8,site);
    xlabel('x0'); ylabel('xfit');
    axis equal
    
    figure
    scatter(meany0,yfit,8,site);
    xlabel('y0'); ylabel('yfit');
    axis equal
    
    figure
    plot(x0,mv_x0,'o');
    xlabel('wn x0'); ylabel('mv x0')
    hold on; axis equal
    plot([0 100],[0 100])
    
    
    figure
    plot(y0,mv_y0,'o')
    xlabel('wn y0'); ylabel('mv y0')
    hold on; axis equal
    plot([0 100],[0 100])
    
    
    figure
    plot(x0,fl_x0,'o');
    xlabel('wn x0'); ylabel('fl x0')
    hold on; axis equal
    plot([0 100],[0 100])
    
    figure
    plot(y0,fl_y0,'o')
    xlabel('wn y0'); ylabel('fl y0')
    hold on; axis equal
    plot([0 100],[0 100])
    
end


% figure
% hold on
% for siteN = 1:32
%     pts = find(site==siteN);
%     plot(histoy(pts),meany0(pts)')
%     if sum(~isnan(meany0(pts)))>1
%         yc = regress(meany0(pts),[ones(size(meany0(pts))) histoy(pts)']);
%         scale(siteN) = yc(2);
%     end
% end

figure
plot(histox,histoy,'o')
axis equal
n_obs=0;
clear obs_name obs

for setup_obs=1:1
    n_obs=n_obs+1; obs_name{n_obs} = 'histox'; obs(:,n_obs) = histox;
    n_obs=n_obs+1; obs_name{n_obs} = 'histoy'; obs(:,n_obs)=histoy;
    n_obs=n_obs+1; obs_name{n_obs} = 'hist sections'; obs(:,n_obs)=histSection;
    n_obs=n_obs+1; obs_name{n_obs} = 'x0'; obs(:,n_obs) = x0;
    n_obs=n_obs+1; obs_name{n_obs} = 'y0'; obs(:,n_obs)=y0;
    n_obs=n_obs+1; obs_name{n_obs} = 'eye'; obs(:,n_obs) = dom_eye;
    n_obs=n_obs+1; obs_name{n_obs} = 'manual type'; obs(:,n_obs) = manual_type;
    n_obs=n_obs+1; obs_name{n_obs} = 'sta amp'; obs(:,n_obs)=A;
    n_obs=n_obs+1; obs_name{n_obs} = 'flash on/off'; obs(:,n_obs) = fl_onoff;
    
    n_obs=n_obs+1; obs_name{n_obs} = 'tminmax'; obs(:,n_obs) = tratio;
    n_obs=n_obs+1; obs_name{n_obs} = 'wx'; obs(:,n_obs) = wx;
    n_obs=n_obs+1; obs_name{n_obs} = 'wy'; obs(:,n_obs)=wy;
    n_obs=n_obs+1; obs_name{n_obs} = 'sta latency'; obs(:,n_obs) = wn_lat;
    n_obs=n_obs+1; obs_name{n_obs} = 'flash latency'; obs(:,n_obs) = fl_lat;
    
    n_obs=n_obs+1; obs_name{n_obs} = 'peak size'; obs(:,n_obs) = fl_sz;
    n_obs=n_obs+1; obs_name{n_obs} = 'fl_lag'; obs(:,n_obs) = fl_lag;
    n_obs=n_obs+1; obs_name{n_obs} = 'offset-onset'; obs(:,n_obs)=-fl_type;
    n_obs=n_obs+1; obs_name{n_obs} = 'sustain'; obs(:,n_obs) = sustain;
    n_obs=n_obs+1; obs_name{n_obs} = 'onfoff corr'; obs(:,n_obs) = fl_onoffcorr;
    
    n_obs=n_obs+1; obs_name{n_obs} = 'mv osi'; obs(:,n_obs) =mv_osi;
    n_obs=n_obs+1; obs_name{n_obs} = 'mv dsi'; obs(:,n_obs) = mv_dsi;
    n_obs=n_obs+1; obs_name{n_obs} = 'peak speed'; obs(:,n_obs) = mv_spd;
    n_obs=n_obs+1; obs_name{n_obs} = 'mv sp ratio'; obs(:,n_obs) = mv_spd_ratio;
    
    n_obs=n_obs+1; obs_name{n_obs} = 'TF ratio'; obs(:,n_obs) = tf_ratio(:,1)';
    n_obs=n_obs+1; obs_name{n_obs} = 'peakSF'; obs(:,n_obs) = peak_sf;
    %n_obs=n_obs+1; obs_name{n_obs} = 'drift amp'; obs(:,n_obs) = sf_amp;
    
    d = driftDSI(:,1)';
    d(d>0.3)=d(d>0.3)+0.1;
    n_obs=n_obs+1; obs_name{n_obs} = 'DSI'; obs(:,n_obs) = d;
    %   n_obs=n_obs+1; obs_name{n_obs} = 'driftDSI1'; obs(:,n_obs) = driftDSI(:,2)';
    
    d = driftOSI(:,1)';
    d(d>0.33)=d(d>0.33)+0.1;
    n_obs=n_obs+1; obs_name{n_obs} = 'OSI'; obs(:,n_obs) = d;
    %  n_obs=n_obs+1; obs_name{n_obs} = 'driftOSI1'; obs(:,n_obs) = driftOSI(:,2)';
    n_obs=n_obs+1; obs_name{n_obs} = 'F1/F0'; obs(:,n_obs) = driftF1F0;
    
    n_obs=n_obs+1; obs_name{n_obs} = 'contrast resp'; obs(:,n_obs) = wn_cr_dom;
    n_obs=n_obs+1; obs_name{n_obs} = 'contrast resp_norm'; obs(:,n_obs)=wn_cr_norm;
    n_obs=n_obs+1; obs_name{n_obs} = 'sbc'; obs(:,n_obs)=(wn_cr_dom<-1*10^3);
    n_obs=n_obs+1; obs_name{n_obs} = 'wn adapt'; obs(:,n_obs)=wn_adapt;
    n_obs=n_obs+1; obs_name{n_obs} = 'spont'; obs(:,n_obs) = dr_spont;
    
end


figure
for i = 1:n_obs;
    subplot(5,7,i)
    [h b]=hist(obs(:,i));
    bar(b,h);
end


%usedcells= find(inside & ~(isnan(obs(:,4)) & isnan(obs(:,16))&isnan(obs(:,24))));
resp =~isnan(obs(:,find(strcmp(obs_name,'sta amp'))))+  ~isnan(obs(:,find(strcmp(obs_name,'offset-onset')))) +...
    ~isnan(obs(:,find(strcmp(obs_name,'peak speed'))))+  ~isnan(obs(:,find(strcmp(obs_name,'F1/F0')))) ;

sbc = wn_cr_dom<-1*10^3;
sbc=sbc';

usedcells = find(manual_type~=0)
%usedcells = find(inside & resp>=1| sbc);
%usedcells = 1:length(sbc);
%usedcells=find(manual_type>=1 & manual_type<=2)
length(usedcells)

obs_used = obs(usedcells,:);
obs_used=replaceNan(obs_used,obs_name,'DSI',0);
obs_used=replaceNan(obs_used,obs_name,'OSI',0);
obs_used=replaceNan(obs_used,obs_name,'mv osi',0);
obs_used=replaceNan(obs_used,obs_name,'mv dsi',0);
obs_used=replaceNan(obs_used,obs_name,'sta amp',0);
obs_used=replaceNan(obs_used,obs_name,'peak speed',3);
obs_used=replaceNan(obs_used,obs_name,'mv sp ratio',0);

obs_norm = (obs_used-repmat(nanmean(obs_used),size(obs_used,1),1))./repmat(nanstd(obs_used),size(obs_used,1),1);
obs_norm(obs_norm>2)=2;
obs_norm(obs_norm<-2)=-2;

figure
for i = 1:n_obs;
    subplot(5,7,i)
    [h b]=hist(obs_norm(:,i));
    bar(b,h);
    set(gca,'Ytick',[])
    title(obs_name{i})
end

labels  ={};
for i = 1:length(obs_name);
    labels{i} = sprintf('%s  %d',obs_name{i},i);
end

figure
pcorr = corrcoef(obs_norm,'rows','pairwise');
pcorr(isnan(pcorr))=0;
imagesc(pcorr)
set(gca,'YTickLabel',labels)
set(gca,'YTick',1:length(obs_name))
title('pairwise')

labels  ={};
for i = 1:length(obs_name);
    labels{i} = sprintf('%s  %d',obs_name{i},i);
end

start_obs=8;
used_obs=start_obs:n_obs;
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'fl_lag'))); %%% fl_lag (1,2 are the only meaningful)
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'contrast resp_norm'))); %%% cr_norm (captured by cr)
%used_obs=used_obs(used_obs~=find(strcmp(obs_name,'driftDSI1'))); %%% drift DSI F1
%used_obs=used_obs(used_obs~=find(strcmp(obs_name,'driftOSI1')));  %%% drift OSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wx'))); %%% drift DSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wy')));  %%% drift OSI F1
%used_obs=used_obs(used_obs~=find(strcmp(obs_name,'sustain')));  %%% sustained / trans from flash; redundant with lag
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'tminmax')));  %%% drift OSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wn adapt')));  %%% drift OSI F1


used_obs=used_obs(used_obs~=find(strcmp(obs_name,'mv sp ratio')));  %%% redundant with speed
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'mv dsi')));  %%% redundant with mv osi
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'mv osi')));  %%% redundant with mv osi
obs_id = find(strcmp(obs_name,'sbc'))
obs_norm(:,obs_id) = obs_norm(:,obs_id)*1.5;
obs_id = find(strcmp(obs_name,'OSI'))
obs_norm(:,obs_id) = obs_norm(:,obs_id)*1.5;
obs_id = find(strcmp(obs_name,'DSI'))
obs_norm(:,obs_id) = obs_norm(:,obs_id)*1.5;
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'sbc')));  %%% redundant with speed


%used_obs=used_obs(used_obs~=find(strcmp(obs_name,'dr spont'))); %%% drift DSI F1

%
% e = eig(pcorr(used_obs,used_obs));
% figure
% plot(e(end:-1:1));
% [v d] = eig(pcorr(used_obs,used_obs));
% figure
% imagesc(v(:,size(v,2):-1:1));
% set(gca,'YTickLabel',labels(used_obs))
% set(gca,'YTick',1:length(obs_name))
% title('pairwise corr')
close all

%%% perform clustering!!!
for rep=2:2;
    for c_rep=6:6
        npca=2*rep+4;
        %npca =6;
        n_clust=c_rep;
        % n_clust=c_rep+5;
        
        phi=1.05                       % fuzzy exponent
        
        maxiter=200;
        toldif=0.0000001;
        
        ndata = length(obs_norm);         % number of data
        % initialise random
        Uinit= initmember(0.1,n_clust,ndata);
        tic
        %data=pc(:,1:npca);
        data = obs_norm(:,used_obs);
        data(isnan(data))=0;
        
        figure
        imagesc(data')
        
        [U, centroid, dist, W, obj] = run_fuzme(n_clust,data,phi,maxiter,1,toldif,0.1,200);
        toc
        %T= kmeans(pc(:,1:7),n_clust,'Replicates',100);
        %T= clusterdata(pc(:,1:6),'maxclust',n_clust,'linkage','ward');
        
        [conf T]= max(U');
        
        for dt=1:2
            [DB(rep,c_rep,dt),CH(rep,c_rep,dt),Dunn(rep,c_rep,dt),KL(rep,c_rep,dt),Han(rep,c_rep,dt),st] = valid_internal_deviation(data,T,dt)
        end
        
        obs_sort = zeros(size(obs_norm));
        obs_sort_mean = zeros(size(obs_norm));
        cell_list = zeros(length(obs_norm),1);
        usort = zeros(size(U));
        n=0;
        for i = 1:max(T);
            members = find(T==i);
            obs_sort(n+1:n+length(members),:) = obs_norm(members,:);
            obs_sort_mean(n+1:n+length(members),:) = repmat(nanmean(obs_norm(members,:),1),length(members),1);
            cell_list(n+1:n+length(members))=members;
            
            usort(n+1:n+length(members),:)=U(members,:);
            n=n+length(members);
        end
        obs_sort(isnan(obs_sort))=-3;
        labels  ={};
        for i = 1:length(obs_name);
            labels{i} = sprintf('%s  %d',obs_name{i},i);
        end
        figure
        imagesc(obs_sort',[-3 3]);
        set(gca,'YTickLabel',labels(1:n_obs))
        set(gca,'YTick',1:length(obs_name))
        figure
        imagesc(obs_sort_mean',[-1.5 1.5]);
        set(gca,'YTickLabel',labels(1:n_obs))
        set(gca,'YTick',1:length(obs_name))
        % figure
        % imagesc(usort')
        
        figure
        imagesc(corrcoef(data(cell_list,:)'),[-1 1]);
        title(sprintf('nclust %d npca %d phi %0.2f',n_clust,npca,phi))
    end
end

figure
imagesc(Dunn(:,:,1));

types = nan(1,length(drift_all));
types(usedcells)=T;

figure
scatter(obs_norm(:,1),obs_norm(:,2),[],T)

[labels cmap clim]= makeColors(types,nan,'qual','Set1');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
colormap(cmap); colorbar; set(gca,'Clim',clim); title(sprintf('type %d',i));


for full_analysis=1:0
dsID = 4; ds = find(types==dsID);
onsustID =2; onsust = find(types == onsustID);
offsustID = 3; offsust= find(types == offsustID);
offtransID =5; offtrans = find(types==offtransID);
sbcID = 1; sbc = find(types==sbcID);
miscwID = 6; miscw = find(types==miscwID);
plist{1} = onsust; plist{2} = offsust; plist{3} = offtrans; plist{4} = ds; plist{5}=miscw; plist{6} = sbc;

figure
for i = 1:length(plist);
    subplot(2,3,i)
    hist(manual_type(plist{i}),1:7)
end

%for i = 1:max(T);
for i = [onsustID offsustID offtransID  dsID   sbcID miscwID];
    id = zeros(size(types));
    id(types==i)=1;
    
    [labels cmap clim]= makeColors(id,nan,'qual','Set1');
    f = plotSections(sections,anatomy,histox(~manual_outside),histoy(~manual_outside),histSection(~manual_outside),labels(~manual_outside,:));
    colormap(cmap); colorbar; set(gca,'Clim',clim); title(sprintf('type %d',i));
    
    f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),id(~manual_outside),[],cbrewer('seq','Greens',64),[0 0.4],100);
    %title(sprintf('type %d',i));
    %subplot(2,3,1); colorbar
end

id = zeros(size(types));
id(types==onsustID | types == offsustID | types ==offtransID)=1;
f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),id(~manual_outside),[],cbrewer('seq','Greens',64),[0 0.7],100);
%title(sprintf('type %d',i));
%subplot(2,3,1); colorbar



sects=[2 3 5];
for section=1:3
    sectionUnits = find(histSection==sects(section) & ~manual_outside');
    Ntotal(section)=length(sectionUnits);
    for type = 1:6;
        N(section,type)=length(intersect(sectionUnits,plist{type}));
        [phat pci] = binofit(N(section,type),Ntotal(section),0.1);
        binoP(section,type)=phat;
        binoCI(section,type,:) = pci;
        err(section,type) = sqrt(Ntotal(section)*phat*(1-phat))/Ntotal(section);
        fraction(section,type) = N(section,type)/Ntotal(section);
    end
end

for tp = 1:6;
    figure
    bar(binoP(:,tp));
    hold on
    %errorbar(1:3,binoP(:,tp),binoP(:,tp)-binoCI(:,tp,1),binoCI(:,tp,2)-binoP(:,tp),'k.','LineWidth',2);
    errorbar(1:3,binoP(:,tp),err(:,tp),'k.','LineWidth',2)
    set(gca,'XTicklabel',{'anterior','mid','posterior'})
    ylabel('fraction'); set(gca,'FontSize',16); set(get(gca,'YLabel'),'Fontsize',16);
end


f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),x0(~manual_outside)',[],jet,[0 90],100);
%subplot(2,3,1)
plot([100 300],[900 900],'Linewidth',4,'Color',[0.5 0.5 0.5])
% subplot(2,3,4); colorbar;
f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),y0(~manual_outside)',[],jet,[-20 40],100);
%subplot(2,3,1)
plot([100 300],[900 900],'Linewidth',4,'Color',[0.5 0.5 0.5])
% subplot(2,3,4); colorbar;



%%% separate ds os

id = zeros(size(types));
id(types' ==dsID & driftDSI(:,1)>0.3)=1;
f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),id(~manual_outside),[],cbrewer('seq','Greens',64),[0 0.3],100);

id = zeros(size(types));
id(types' ==dsID & driftDSI(:,1)<0.3)=1;
f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),id(~manual_outside),[],cbrewer('seq','Greens',64),[0 0.4],100);

[labels cmap clim]= makeColors(id,nan,'qual','Set1');
f = plotSections(sections,anatomy,histox(~manual_outside),histoy(~manual_outside),histSection(~manual_outside),labels(~manual_outside,:));
colormap(cmap); colorbar; set(gca,'Clim',clim); title(sprintf('type %d',i));


ncells = [length(offtrans) length(onsust) length(offsust) length(ds) length(sbc) length(miscw) sum(~manual_outside' & isnan(types))];
figure
bar(ncells/sum(ncells));
grouplabels = {'tOff','sOn','sOff','DS/OS','supp','slow','unresp'};
set(gca,'XTicklabel',grouplabels)
rotateticklabel(gca,90)
ylabel('fraction')
set(get(gca,'Ylabel'),'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
ylim([0 0.28])

figure
for i = 1:length(plist);
    subplot(2,3,i)
    hist(manual_type(plist{i}),1:7)
end

%%% Figures for lgn paper

%%% figure 1 - ephys & sections

%%% correlation
obs_sort = zeros(size(data));
obs_sort_mean = zeros(size(data));
cell_list = zeros(length(obs_norm),1);
usort = zeros(size(U));
n=0;
for i = [onsustID offsustID offtransID  dsID   sbcID miscwID];
    members = find(T==i);
    obs_sort(n+1:n+length(members),:) = obs_norm(members,used_obs);
    obs_sort_mean(n+1:n+length(members),:) = repmat(nanmean(obs_norm(members,used_obs),1),length(members),1);
    cell_list(n+1:n+length(members))=members;
    
    usort(n+1:n+length(members),:)=U(members,:);
    n=n+length(members);
end
%obs_sort(isnan(obs_sort))=-3;
labels  ={};
for i = 1:length(used_obs);
    labels{i} = sprintf('%s',obs_name{used_obs(i)});
end
figure
imagesc(obs_sort',[-3 3]);
set(gca,'YTickLabel',labels(1:length(used_obs)))
set(gca,'YTick',1:length(obs_name))
colormap(flipud(cbrewer('div','RdBu',64)));

figure
imagesc(obs_sort_mean',[-1 1]);
set(gca,'YTickLabel',labels(1:length(used_obs)))
set(gca,'YTick',1:length(obs_name))
%colormap(gray)
colormap(flipud(cbrewer('div','RdBu',64)));
set(gca,'Xtick',50:50:200); xlabel('cell #');
colorbar
set(get(gca,'Xlabel'),'FontSize',16); set(gca,'Fontsize',14);

%%% correlation
figure
imagesc(corrcoef(data(cell_list,:)'),[-1 1]);
%colormap gray
colormap(flipud(cbrewer('div','RdBu',64)));
colorbar
set(gca,'Xtick',50:50:200);
set(gca,'Ytick',50:50:200);
set(gca,'Fontsize',14)

[c p]=corrcoef(data(cell_list,:)');
figure
imagesc(p<0.05);
mean(c(163:176,163:176));

c= corrcoef(data');
order= [offtransID onsustID offsustID dsID   sbcID miscwID];
for i = 1:6
    groupc = c(find(T==order(i)),find(T==order(i)));
    cc(i,1) = mean(groupc(:));
    ccstd(i,1) = std(groupc(:))/sqrt(sqrt(0.5*length(groupc(:))));
    
    groupc = c(find(T==order(i)),find(T~=order(i)));
    cc(i,2) = mean(groupc(:));
    ccstd(i,2) = std(groupc(:))/sqrt(sqrt(0.5*length(groupc(:))));
end

ticksize=14; labelsize=20; stdwidth=2;

grouplabels = {'tOff','sOn','sOff','DS/OS','supp','slow'};

figure
bar(cc(:,1));
hold on
bar(cc(:,2),'r');
errorbar(cc(:,2),ccstd(:,2),'k.','LineWidth',2);
errorbar(cc(:,1),ccstd(:,1),'k.','LineWidth',2);
axis([0 7 -0.25 0.75])
ylabel('correlation');
set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
legend('within group', 'btw groups')
set(gca,'Xticklabel',grouplabels);
set(get(gca,'Xlabel'),'FontSize',16);

%%% size suppresion
figure
fl_cc = fl_supp([plist{1},plist{2},plist{3}]);
fl_cc(fl_cc<0)=0;
hist(1-fl_cc,0:0.2:1)
set(gca,'Ytick',13.5:13.5:45); set(gca,'Yticklabel',0.1:0.1:0.3);
xlabel('size supp');
ylabel('fraction')
set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
xlim([-0.2 1.2])

figure
bandpass = lowpass([plist{1},plist{2},plist{3}])./drift_amp([plist{1},plist{2},plist{3}]);
bandpass(bandpass>1)=1;
hist(1-bandpass,0:0.2:1)
set(gca,'Ytick',13.5:13.5:45); set(gca,'Yticklabel',0.1:0.1:0.3);
xlabel('bandpass index');
ylabel('fraction')
set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
xlim([-0.2 1.2])



%%% bar graphs
figure
subplot(3,3,1)
barparams(driftOSI,plist);
ylabel('OSI','FontSize',10);
subplotAxes(gca)


subplot(3,3,2);
hold off
barparams(A*8,plist);
ylabel('sta amp', 'Fontsize',10);
ylim([-1 1])
subplotAxes(gca)
% subplot(3,3,2);
% barparams(sustain,plist);
% ylabel('sust','FontSize',10);
% subplotAxes(gca)

subplot(3,3,3);
barparams_med(driftF1F0*2,plist);
ylabel('F1/F0','FontSize',10);
subplotAxes(gca)


subplot(3,3,5);
speeds = [10 20 40 80 160 nan];
mv_spd(isnan(mv_spd))=6;
mv_spd_fix = speeds(mv_spd);
barparams(mv_spd_fix,plist);
ylabel('peak speed deg/sec','Fontsize',10)
subplotAxes(gca)


subplot(3,3,9);
sfs = [0 0.01 0.02 0.04 0.08 0.16 0.32 nan];
sf_inds = peak_sf;
sf_inds(peak_sf ==0) = 7;
sf_inds(isnan(peak_sf))=8;
barparams(sfs(sf_inds),plist);
ylabel('peak SF (cpd)','FontSize',10)
ylim([0 0.14])
subplotAxes(gca)

% subplot(3,3,7);
% barparams(tf_ratio(:,1),plist);
% ylabel('TF ratio', 'Fontsize',10);

subplot(3,3,6);hold off
barparams_med(wn_cr_dom/600,plist,6);
ylabel('wn resp sp/sec');
ylim([-8 8])
subplotAxes(gca,6)


subplot(3,3,7);hold off
barparams_med(dr_spont,plist,6);
ylabel('spont (sp/sec)', 'Fontsize',10);
ylim([0 11])
subplotAxes(gca,6)


subplot(3,3,8);hold off
barparams(fl_onoffcorr,plist);
ylabel('on/off corr', 'Fontsize',10);
subplotAxes(gca)

subplot(3,3,8);hold off
barparams(tratio,plist);
ylabel('sta transience', 'Fontsize',10);
subplotAxes(gca)



% subplot(3,3,8);
% barparams(wx,plist);
% ylabel('wx', 'Fontsize',10);
%
% subplot(3,3,9);
% barparams(wy,plist);
% ylabel('wy', 'Fontsize',10);


%democells = [57 89 199]; %sbc
%democells = 89
% democells=199
%  democells = [124 129 157 182 185 231] %% ds os
% democells = [124 157 182]
%  democells = [ 156] %sOff
% democells = [108 ] %sOff
% democells = [119 ]; %sOn
democells = [181 236]; %tOff
% democells = [181]; %tOff
xedge = 64-32;
ticksize=14; labelsize=20; stdwidth=2;
for c = 1:length(democells);
    figure
    i=democells(c); cell_n=i;
    sta = squeeze(wn_all(i,dom_eye(i)).svd_xy(1,:,:));
    imagesc(sta,1.2*[-max(abs(sta(:))) max(abs(sta(:)))]);
    colormap(gray(256));
    xlim([xedge+1 xedge+64]);
    ylim([1 64])
    %       set(gca,'Xtick',xedge:10/wn_degperpix:xedge+64);
    %       set(gca,'Xticklabel',{'0','10','20','30','40','50'});
    
    hold on
    plot([35 35+10/wn_degperpix],[5 5],'w','Linewidth',12)
    %        set(gca,'Ytick',1:10/wn_degperpix:64);
    %       set(gca,'Yticklabel',{'0','10','20','30','40','50'});
    %       xlabel('deg');
    axis off
    axis square
    
    
    
    dt=0.05;
    tmax=0.25;
    figure
    stat= squeeze(wn_all(i,dom_eye(i)).sta_t);
    
    if ~isnan(stat)
        if abs(min(sta(:)))>max(sta(:))
            stat = -stat;
        end
        plot(stat,'Linewidth',2)
        ylim(1.2*[-max(abs(stat)) max(abs(stat))]);
        xlim([1 16])
        set(gca,'Xtick',(0:dt:tmax)*60+1);
        set(gca,'Xticklabel',{'0', '50','100','150','200','250'});
        set(gca,'Ytick',[-0.5 -0.25 0 0.25 0.5])
        xlabel('msec')
        set(get(gca,'Xlabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    end
    
    figure
    f= fl_all(i);
    for r = 1:2
        hold on
        if r==1
            plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'Linewidth',2)
        else
            plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'r','Linewidth',2)
        end
    end
    yl = get(gca,'YLim');
    plot([0.25 0.5],0.9*[yl(2) yl(2)],'g','Linewidth',4);
    xlim([0.15 f.onset_bins(end-1)]);
    set(gca,'Xtick',0.15:0.1:f.onset_bins(end-1));
    set(gca,'Xticklabel',((0.15:0.1:f.onset_bins(end-1))-0.15)*1000)
    xlabel('msec')
    ylabel('sp/sec');
    set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    
    figure
    plot(mv_sp_tune(cell_n,:),'Linewidth',2);
    if ~isnan(mv_sp_tune(cell_n))
        axis([1 5 min(0,min(mv_sp_tune(cell_n,:))) 1.2*max(mv_sp_tune(cell_n,:))])
    end
    set(gca,'XTick',1:5); set(gca,'XtickLabel',[10 20 40 80 160]); xlabel('deg/sec');
    ylabel('sp/sec')
    set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    
    d = drift_all(cell_n);
    color = {'b','r'};
    figure
    for r=1:2  %%% tempfreq
        f=1;
        drift_os = squeeze(d.orient_tune(r,f,:));
        h=polar((0:8)*pi/4,[drift_os' drift_os(1)],color{r});
        set(h,'Linewidth',stdwidth)
        hold on
    end
    
    
    figure
    for r=1:2
        drift_sf= squeeze(d.sf_tune(r,f,:));
        
        hold on
        plot(drift_sf,color{r},'Linewidth',2);
        axis([1 7 min(0,min(drift_sf)) drift_amp(cell_n)+1])
        set(gca,'XTick',1:7); set(gca,'Xticklabel',[0 0.01 0.02 0.04 0.08 0.16 0.32]);
        xlabel('SF (cpd)'); ylabel('sp/sec')
    end
    set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    legend('2Hz','8Hz')
    
    figure
    plot(fl_sztune(i,:),'Linewidth',stdwidth);
    xlim([1 6]); ylim([0 1.2*max(fl_sztune(i,:))]);
    set(gca,'Xticklabel',{'2','4','8','16','32','full'}); xlabel('spot size'); ylabel('sp/sec');
    set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    
    figure
    plot((1-cos((1:20)*2*pi/20))*0.5,squeeze(wn_all(i,1).crf)/600,'Linewidth',stdwidth);
    xlabel('contrast'); ylabel('sp/sec');
    set(gca,'Xtick',[0 0.25 0.5 0.75 1]); ylim([0 1.2*max(wn_all(i,1).crf/600)]);
    set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);
    
    
end

figure
plot(dr_spont(~manual_type==0), wn_cr_dom(~manual_type==0)/600,'o')
hold on
%plot(dr_spont(plist{6}), wn_cr(plist{6})/600,'ro')
plot(dr_spont(manual_type==4), wn_cr_dom(manual_type==4)/600,'go')

figure
plot(wx([plist{1} plist{2} plist{3}]),wy([plist{1} plist{2} plist{3}]),'o')

labelsize=20; ticksize=16;
figure
hist(1.17*wx([plist{1} plist{2} plist{3}]),1:1:14)
xlabel('radius (deg)');
ylabel('fraction');
npts = length([plist{1} plist{2} plist{3}])
set(gca,'Xtick',2:2:14)
set(gca,'Ytick',(0:0.1:0.4)*npts); set(gca,'YTickLabel',0:0.1:0.4);
axis([0 15.5 0 0.35*npts])
set(get(gca,'Xlabel'),'FontSize',labelsize); set(get(gca,'Ylabel'),'FontSize',labelsize); set(gca,'Fontsize',ticksize);


figure
os = find(types==dsID);
plot(driftOSI(os,1),driftDSI(os,1),'o')
xlabel('OSI'); ylabel('DSI');
axis equal
axis square


id = zeros(1,length(types));
id(types' ==dsID & (sin(driftOSItheta(:,1)).^2)>0.5)=1;
f = plotFiltSections(sections,anatomy,round(histox(~manual_outside)),round(histoy(~manual_outside)),histSection(~manual_outside),id(~manual_outside),[],cbrewer('seq','Greens',64),[0 0.4],100);

hrange = 0.05:0.1:1
figure
hist(driftOSI(~isnan(types) ),hrange);
hold on
hist(driftOSI(types==dsID),hrange,'r');

driftDSI(driftDSI(:,1)>0.4 & types'~=dsID)=nan;
hrange = 0.05:0.1:1
figure
hist(driftDSI(~isnan(types) ),hrange);
hold on
hist(driftDSI(types==dsID),hrange,'r');
sum(hist(driftDSI(~isnan(types) ),hrange))


grouplabels = {'sOn','sOff','tOff','DS/OS','slow','supp'};

datatable=[];
for i = 1:6
    datatable{1,i+1}=grouplabels{i};
end
for i = 1:size(obs,2);
    datatable{i+1,1}=obs_name{i};
    for c = 1:6
        if i==8
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',9*nanmean(obs(plist{c},i)),9*nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
        elseif i==13
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',16.6*nanmean(obs(plist{c},i)),16.6*nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
        elseif i==14
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',1000*nanmean(obs(plist{c},i))/20,(1000/20)*nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
        elseif i==15
            szs= [2 4 8 16 32 64 nan];
            sz = obs(:,i); sz(isnan(sz))=7;
            sz_all = szs(sz);
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',nanmean(sz_all(plist{c})),nanstd(sz_all(plist{c}))/sqrt(length(plist{c})));
        elseif i==22
            speeds = [10 20 40 80 160 nan];
            mv_spd(isnan(mv_spd))=6;
            mv_spd_fix = speeds(mv_spd);
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',nanmean(mv_spd_fix(plist{c})),nanstd(mv_spd_fix(plist{c}))/sqrt(length(plist{c})));
        elseif i == 25
            sfs = [0 0.01 0.02 0.04 0.08 0.16 0.32 nan];
            sf_inds = peak_sf;
            sf_inds(peak_sf ==0) = 7;
            sf_inds(isnan(peak_sf))=8;
            sfvalue = sfs(sf_inds);
            datatable{i+1,c+1} = sprintf('%0.2f +/- %0.3f',nanmean(sfvalue(plist{c})),nanstd(sfvalue(plist{c}))/sqrt(length(plist{c})));
        elseif i==28
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',2*nanmedian(obs(plist{c},i)),2*nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
        elseif i==29
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',nanmedian(obs(plist{c},i))/600,(1/600)*nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
            
        else
            datatable{i+1,c+1} = sprintf('%0.1f +/- %0.2f',nanmedian(obs(plist{c},i)),nanstd(obs(plist{c},i))/sqrt(length(plist{c})));
        end
    end
end
xlswrite('datatableLGN2.xls',datatable)


figure
histrange = 0:45:330
rose([ driftOSItheta(os,1) ;driftOSItheta(os,1)+pi],histrange*pi/180)
set(gca,'LineWidth',4)
set(gca,'Fontsize',16)

figure
histrange = 0:45:180
h=hist(round(driftOSItheta(os,1)*180/pi),histrange)
hfix= h;
hfix(1)=h(1)+h(end);
hfix(end) = h(1)+h(end);
polar(histrange*pi/180,hfix)
hold on
polar(histrange*pi/180+pi,hfix);

ds = find(types' == dsID & driftDSI(:,1)>0.3)
figure
histrange = 0:45:330;
rose(round(driftDSItheta(ds,1)),histrange*pi/180)

h=hist(round(driftDSItheta(ds,1)*180/pi),histrange);
hfix= h;
hfix(1)=h(1)+h(end);
hfix(end) = h(1)+h(end);
polar(histrange*pi/180,hfix)

c = find(~isnan(types));
normY= nan(size(histox));
depthY= nan(size(histox));
for cell_n = c;
    
    paxxy = sections(histSection(cell_n)).coords;
    tracexy = anatomy(site(cell_n)).LGN;  %%% coords are reversed
    widthX = max(tracexy(:,2))-min(tracexy(:,2));
    %normX(cell_n) = (histox(cell_n)-min(tracexy(:,2))) / widthX;
    ypts = find(abs((tracexy(:,2))-(histox(cell_n)))<10);
    if length(ypts)>1
        normY(cell_n) = (histoy(cell_n)-min(tracexy(ypts,1))) / (max(tracexy(ypts,1)) - min(tracexy(ypts,1)));
        depthY(cell_n) = max(tracexy(ypts,1)) -histoy(cell_n);
    else
        normY(cell_n)=-50;
        depthY(cell_n)=-50;
    end
end
figure
allY=hist(depthY(~manual_outside & histSection'~=5),[150 450])
bar(allY)
hold on
dsY=hist(depthY( find(types==dsID & histSection~=5)),[150 450]);
bar(dsY)
figure
bar(dsY./allY)

range = [0:0.25:1]
figure
allY=hist(normY,range);
hist(normY,range);
hold on
dsY=hist(normY( find(types==dsID)),range);
hist(normY( find(types==dsID)),range);
figure
bar(dsY./allY)
end %%% full analysis

psfilename = 'c:/data/allLGN.ps'
if exist(psfilename,'file')==2;delete(psfilename);end
figure
imagesc(obs_sort_mean');
labels  ={};
for i = 1:length(obs_name);
    labels{i} = sprintf('%s  %d',obs_name{i},i);
end
set(gca,'YTickLabel',labels(1:n_obs))
set(gca,'YTick',1:length(obs_name))
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

%for c = 1:length(cell_list);
for c = 1:length(wn_all)
    
    
    %cell_n = usedcells(cell_list(c));
    cell_n=c;
    
    i = cell_n;
    
    figure
    set(gcf,'position',[100 50 850 950])
    subplot(5,4,1);
    
    if ~isempty(wn_all(cell_n,dom_eye(cell_n)).N)
        for j=1:3
            subplot(5,4,j);
            imagesc(squeeze(wn_all(i,dom_eye(cell_n)).svd_xy(j,:,:)));
            axis off
        end
        
        subplot(5,4,5);
        imagesc(squeeze(wn_all(i,dom_eye(cell_n)).fitsta));
        axis off
        title(sprintf('x=%d y=%d wx=%d wy=%d',round(x0(cell_n)),round(y0(cell_n)),round(wx(cell_n)),round(wy(cell_n))));
        subplot(5,4,6);
        plot(squeeze(wn_all(i,dom_eye(cell_n)).sta_t));
        title(sprintf('tminmax = %0.2f',tratio(cell_n)));
        subplot(5,4,7);
        plot(1-cos((1:20)*2*pi/20),squeeze(wn_all(i,1).crf));hold on
        if ~isempty(wn_all(i,2).crf)
            plot(1-cos((1:20)*2*pi/20),squeeze(wn_all(i,2).crf),'r');
        end
        title(sprintf('cnrm=%0.1f ad=%0.1f dom=%d',wn_cr_norm(cell_n),wn_adapt(cell_n),dom_eye(cell_n)))
    end
    subplot(5,4,1);
    if ismember(cell_n,usedcells)
        title(sprintf('cell %d; type %d',c,T(find(usedcells==c))));
    else
        title(sprintf('cell %d; type nan',c))
    end
    subplot(5,4,2);
    title(sprintf('filenum %d',site(cell_n)));
    
    subplot(5,4,3);
    title(sprintf('ch%d cl%d',cell_id(cell_n,1),cell_id(cell_n,2)));
    
    f= fl_all(i);
    subplot(5,4,9);
    if ~isempty(f.sta_diff);
        imagesc(f.sta_diff(:,:,f.lag+1)');
    end
    title(sprintf('flash amp=%0.1f',fl_amp(cell_n)))
    axis off
    
    if ~isempty(f.N)
        subplot(5,4,10);
        for r = 1:2
            hold on
            if r==1
                plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)))
            else
                plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'r')
            end
        end
        yl = get(gca,'YLim');
        plot([0.25 0.5],0.9*[yl(2) yl(2)],'g')
        xlim([0 f.onset_bins(end-1)])
        title(sprintf('tp=%d sust=%0.1f c=%0.1f',fl_type(cell_n),sustain(cell_n),fl_onoffcorr(cell_n)));
        
        subplot(5,4,11);
        plot(fl_sztune(i,:));
        xlim([1 6])
        set(gca,'XTick',[]); xlabel('size');
    end
    
    if isstruct(mv_all(cell_n))
    if ~isempty(mv_all(cell_n).N)
        subplot(5,4,13);
        imagesc(squeeze(mv_all(cell_n).sta_diff(:,:,mv_all(cell_n).lag+1))');
        axis off
        title(sprintf('moving amp=%0.1f',mv_amp(cell_n)))
        
        subplot(5,4,14);
        plot(mv_sz_tune(cell_n,:));
        if ~isnan(mv_sz_tune(cell_n,:))
            axis([1 3 min(0,min(mv_sz_tune(cell_n,:))) 1+max(mv_sz_tune(cell_n,:))])
        end
        set(gca,'XTick',[]); xlabel('size');
        
        subplot(5,4,15);
        plot(mv_sp_tune(cell_n,:));
        if ~isnan(mv_sp_tune(cell_n))
            axis([1 5 min(0,min(mv_sp_tune(cell_n,:))) 1+max(mv_sp_tune(cell_n,:))])
        end
        set(gca,'XTick',[]); xlabel('speed');
        title(sprintf('sp=%0.1f sp_ratio=%0.1f',mv_spd(cell_n),mv_spd_ratio(cell_n)));
        subplot(5,4,16);
        plot(mv_os_tune(cell_n,:));
        title(sprintf('os=%0.1f ds=%0.1f',mv_osi(cell_n),mv_dsi(cell_n)));
        if~isnan(mv_os_tune(cell_n))
            axis([1 8 min(0,min(mv_os_tune(cell_n,:))) 1+max(mv_os_tune(cell_n,:))]);
        end
        set(gca,'XTick',[]); xlabel('theta');
    end
    else
        subplot(5,4,13)
        bar(burst_fraction(cell_n,:));
        ylim([0 0.4]);
        title(sprintf('burst %0.2f %0.2f',burst_fraction(cell_n,1),burst_fraction(cell_n,2)));
        subplot(5,4,14);
        plot(stop_crf(cell_n,:),'r');
        hold on
        plot(move_crf(cell_n,:),'b');
        yl = get(gca,'ylim');
        ylim([0 yl(2)]);
        
        if length(movement_all(cell_n).speed)>1
            subplot(5,4,15)
        plot(movement_all(cell_n).mvlfp_tdata,movement_all(cell_n).speed,'g');
        title(sprintf('moving = %0.0f %%',100*sum(movement_all(cell_n).speed>1) /length(movement_all(cell_n).speed)));
        xlim([0 max(movement_all(cell_n).mvlfp_tdata)])
        
        subplot(5,4,16)
          plot(movement_all(cell_n).freqs, squeeze(movement_all(cell_n).mv_lfp));
        ylim([0 1.5*prctile(movement_all(cell_n).mv_lfp(1,:),95)])
        xlim([0 90])
        end
        
        %%% evoked ratio?
    end
    if ~isnan(histox(cell_n))
    subplot(5,4,4)
    plot(sections(histSection(cell_n)).coords(:,1),sections(histSection(cell_n)).coords(:,2),'b.'); hold on
    plot(originalX(cell_n),originalY(cell_n),'r*');
    hold on
    plot(alignedX(cell_n),alignedY(cell_n),'g*')
    title(sprintf('Pax AP = %d',histSection(cell_n)));
    axis([-500 500 -500 500]);
    axis off
    
    subplot(5,4,8);
    lgnxy = anatomy(site(cell_n)).LGN;
    plot(lgnxy(:,2),lgnxy(:,1),'b.'); hold on
    plot(originalX(cell_n),originalY(cell_n),'r*');
    etrode = anatomy(site(cell_n)).Electrode;
    plot(etrode(:,2),etrode(:,1),'k','LineWidth',1);
    axis([-500 500 -500 500]);
    axis off
    title(sprintf('monitor %0.0f',offsetX(cell_n)));
    end
    
    subplot(5,4,12);
    plot(wvform(cell_n,:));
    set(gca,'Xtick',[]);
    xlabel('wvform');
    
    
    color = {'b','r'};
    d = drift_all(cell_n);
    for r=1:2  %%% tempfreq
        for f=1:2  %%% F0 F1
            drift_os = squeeze(d.orient_tune(r,f,:));
            
            subplot(5,4,16+f)
            hold on
            plot(drift_os,color{r});
            title(sprintf('os%0.1f %d ds%0.1f %d',driftOSI(cell_n,f),round(driftOSItheta(cell_n,f)*180/pi),driftDSI(cell_n,f), round(driftDSItheta(cell_n,f)*180/pi)));
            axis([1 8 min(0,min(drift_os)) drift_amp(cell_n)+1])
            set(gca,'XTick',[])
            xlabel('theta')
            
            drift_sf= squeeze(d.sf_tune(r,f,:));
            subplot(5,4,18+f)
            hold on
            plot(drift_sf,color{r});
            if f==2
                title(sprintf('peak sf %d',peak_sf(cell_n)));
            else
                title(sprintf('sf amp %0.1f F1F0 %0.1f',drift_amp(cell_n),driftF1F0(cell_n)));
            end
            axis([1 7 min(0,min(drift_sf)) drift_amp(cell_n)+1])
            set(gca,'XTick',[])
            xlabel('SF')
        end
    end
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    close(gcf)
    
end
ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);


