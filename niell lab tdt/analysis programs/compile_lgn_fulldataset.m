close all
clear all
n_obs=0;
afile = {    'C:\data\lgn rf project_new0824\lgn_analysis\030912_rec3_analysis3_final.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\030912_rec4_analysis_4final.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031012_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031012_analysis2check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031512_analysis2check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031512_analysis3check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031612_analysis1checkb.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031612_analysischeck2b.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\031612_analysis3check.mat', ...
    'C:\data\lgn rf project_new0824\lgn_analysis\052512_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\060112_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\060712_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\060812_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\061212_analysis1check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\061212_analysis2check.mat',...
    'C:\data\lgn rf project_new0824\lgn_analysis\071612_analysis1new.mat', ...
    'C:\data\lgn rf project_new0824\lgn_analysis\071712_analysis1new.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\071712_analysis2new.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\071812_analysis1new.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\072012_analysis1newst.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\072012_analysis2new.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\072012_analysis3.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\072412_analysis1.mat' ...
    'C:\data\lgn rf project_new0824\lgn_analysis\072512_analysis1.mat' ...
     'C:\data\lgn rf project_new0824\lgn_analysis\072512_ANALYSIS2.mat' ...
      'C:\data\lgn rf project_new0824\lgn_analysis\080212_analysis1.mat' ...
       'C:\data\lgn rf project_new0824\lgn_analysis\080212_analysis3.mat' ...
        'C:\data\lgn rf project_new0824\lgn_analysis\080312_analysis2.mat' ...
          'C:\data\lgn rf project_new0824\lgn_analysis\080312_analysis3bnew.mat' };


histfile = { 'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn015_P3_A53.tif', ...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn015_P4_A51.tif', ...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn016_P5_A53.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn016_P6_A53.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn017_P4_A53.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn017_P7_A52.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn018_P3_A54.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn018_P4_A51.tiff',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn018_P5_A52.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\DMN020P1F50.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn021P7F50.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn022P10F48.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn023P5F53.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\DMN024P9F53.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\DMN024P10F51.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn026s10p4f53r3.tif', ...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn026s10p5f53r4.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn027S10P1F49R1.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn028S11P5F49R2.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn028S12P3F47R1.tif',...
    'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn029S5P5F55R2.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn029S6P2F54R1tif.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn030S9P2F52R1.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn030S9P3F53R2.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn030S9P4F52R3.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn031S6P3F52R1.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn032S5P6F54R2.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn032S7P5F51R1.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn033S4P3F1R1.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn033S4P5F51R3.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn034S6P3F51R2.tif',...
     'C:\data\lgn rf project_new0824\histology data\ANATOMY2\dmn034S6P4F51R3.tif',...
    };

%PaxPos=[54 53 50 53 53 53 52 50 53 54 50 50 48 53 53 51 53];
PaxPos=[54 53 50 53 53 52 53 53 50 54 50 50 48 53 53 51 53 53 53];
%  AllenPos=[82 81 77 81 81 83 76 81 80 77 81 82 77 77 75 81 81 79 81 80];

n=0;
%PaxPos=[55 53 50 51 51 53 49 48 53 54];
[anatomy sections]=LGN_histo(histfile, PaxPos);

numsites = 32*ones(1,length(histfile));
numsites(12)=16;


for i = 1:length(afile);
    i
    clear displayOffset drift mv fl wn
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
    mv_all(cell_range)=mv;
    fl_all(cell_range)=fl;
    site(cell_range)=i;
    cell_id(cell_range,:) = cells;
    wvform(cell_range,:) = wv';
    if ~exist('displayOffset')
        display(sprintf('need to add display for file %d',i))
        %input('press return to continue')
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
        x(j)=anatomy(i).siteXY(1,peakchan(j));
        y(j)=anatomy(i).siteXY(2,peakchan(j));
        z(j)=anatomy(i).AP;
        s(j)=anatomy(i).section;
        
    end
    
    histox(cell_range)=x;
    histoy(cell_range)=y;
    histoz(cell_range)=z;
    histSection(cell_range)=s;
    
end

  for cell_n=1:length(histox)
  %for cell_n=1:1  
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
    ypts = find(abs((paxxy(:,2))-(alignedX(cell_n)))<10);
    if length(ypts)>1
        alignedY(cell_n) = min(paxxy(ypts,2))+normY(cell_n)*(max(paxxy(ypts,2))-min(paxxy(ypts,2)));
    else
        alignedY(cell_n)=0;
    end
  end
    

inside= ones(n,1);
for cell_n= 1:n
    sec = sections(histSection(cell_n)).coords;
    y= sec(find(round(sec(:,1))==round(histox(cell_n))),2);
    if isempty(y) || histoy(cell_n)>max(y) || histoy(cell_n)<min(y)
        inside(cell_n)=0;
    end
end


labels = zeros(n,3);
labels(:,3)=1;
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);


wn_degperpix = wn(end).degperpix;

clear p

n_sta=0;

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
                
                if z>6
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

ind=0;
tmin=zeros(1,cell_n); tmax = tmin; wn_lat=tmin;

x0 = nan(n,1); y0 = nan(n,1); wx=nan(n,1),wy=nan(n,1);A = nan(n,1); wn_resp = nan(n,1); wn_phase=nan(n,1);
tmin= nan(n,1); tmax=nan(n,1); wn_lat = nan(n,1); wn_cr_norm=nan(n,1);wn_adapt=nan(n,1); dom_eye=ones(n,1)

dom_eye(wn_cr(:,2)>1000 & wn_cr(:,2)>wn_cr(:,1)& wn_cr(:,1)>-10^3)=2;
figure
scatter(wn_cr(:,1),wn_cr(:,2),6,dom_eye)

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
                y0(cell_n)= (wn_all(cell_n,eye).sta_fit(3) - 38)*wn_degperpix -offsetY(cell_n);  %%% top of monitor = -38, bottom = 38,  so moving upward means more negative
                y0(cell_n) = 90 - y0(cell_n);  %%% invert so y goes from bottom up
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
            end
        end
    end
end

n_obs=n_obs+1; obs_name{n_obs} = 'histox'; obs(:,n_obs) = histox;
n_obs=n_obs+1; obs_name{n_obs} = 'histoy'; obs(:,n_obs)=histoy;
n_obs=n_obs+1; obs_name{n_obs} = 'hist sections'; obs(:,n_obs)=histSection;
n_obs=n_obs+1; obs_name{n_obs} = 'x0'; obs(:,n_obs) = x0;
n_obs=n_obs+1; obs_name{n_obs} = 'y0'; obs(:,n_obs)=y0;
n_obs=n_obs+1; obs_name{n_obs} = 'eye'; obs(:,n_obs) = dom_eye;
n_obs=n_obs+1; obs_name{n_obs} = 'sta amp'; obs(:,n_obs)=A;

figure
plot (histox,x0,'o');
xlabel('histology X');
ylabel('RF X');


figure
plot (histoy,y0,'o');
xlabel('histology Y');
ylabel('RF Y');


labels=ones(n,3);
goodfit = find(~isnan(x0));
labels(goodfit,1)=0;
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
title('wn sta rfs');



[labels cmap clim]= makeColors(x0,nan,'div','RdBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
colormap(cmap); colorbar; set(gca,'Clim',clim); title('x0');

[labels]= makeColors(y0,nan,'div','RdBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
colormap(cmap); colorbar; set(gca,'Clim',clim); title('y0');

tratio=-tmin./tmax;
tratio(tratio>1)=1;
labels= makeColors(tratio,nan,'div','RdBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
title('min max ratio');


n_obs=n_obs+1; obs_name{n_obs} = 'tminmax'; obs(:,n_obs) = tratio;
n_obs=n_obs+1; obs_name{n_obs} = 'wx'; obs(:,n_obs) = wx;
n_obs=n_obs+1; obs_name{n_obs} = 'wy'; obs(:,n_obs)=wy;
n_obs=n_obs+1; obs_name{n_obs} = 'wn_latency'; obs(:,n_obs) = wn_lat;

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
        
        fl_lat(cell_n)=onset_lat(cell_n,fl_onoff(cell_n));
        fl_sust(cell_n) =onset_sust(cell_n,fl_onoff(cell_n));
        fl_onset_amp(cell_n) = onset(cell_n,fl_onoff(cell_n));
        
        fl_lag(cell_n) = f.lag;
        fl_onoff(cell_n) = 2*(fl_onoff(cell_n)-1.5);
        if fl_lag(cell_n)==1  %%% if a cell responds with 1 frame delay, it is responding to offset
            fl_onoff(cell_n)=-1*fl_onoff(cell_n);
        end
    end
end

fl_type= nan(size(fl_lag));
fl_type(fl_lag==0)=1;
fl_type(fl_lag==1)=2;
labels= makeColors(fl_type,nan);
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
title('flash type')

%n_obs=n_obs+1; obs_name{n_obs} = 'fl size supp'; obs(:,n_obs) = fl_supp;
fl_thresh=1;

fl_sz(fl_amp<fl_thresh)=nan;
fl_lag(fl_amp<fl_thresh)=nan;
fl_type(fl_amp<fl_thresh)=nan;
fl_sust(fl_amp<fl_thresh)=nan;
fl_onoff(fl_amp<fl_thresh)=nan;
fl_lat(fl_amp<fl_thresh)=nan;
fl_onoffcorr(fl_amp<fl_thresh)=nan;


n_obs=n_obs+1; obs_name{n_obs} = 'fl peak size'; obs(:,n_obs) = fl_sz;
n_obs=n_obs+1; obs_name{n_obs} = 'fl_lag'; obs(:,n_obs) = fl_lag;
n_obs=n_obs+1; obs_name{n_obs} = 'fl_type'; obs(:,n_obs)=fl_type;
n_obs=n_obs+1; obs_name{n_obs} = 'fl_sust'; obs(:,n_obs) = fl_sust;
n_obs=n_obs+1; obs_name{n_obs} = 'fl_on/off'; obs(:,n_obs) = fl_onoff;
n_obs=n_obs+1; obs_name{n_obs} = 'fl latency'; obs(:,n_obs) = fl_lat;
n_obs=n_obs+1; obs_name{n_obs} = 'onfoff corr'; obs(:,n_obs) = fl_onoffcorr;



labels= makeColors(fl_type,nan,'seq','GnBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
title('fl_type')

fl_sust_norm = fl_sust;
fl_sust_norm(fl_sust_norm<0)=0;
fl_sust_norm(fl_onset_amp<2)=nan;
labels= makeColors(fl_sust_norm);
f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
title('normalized transience')

clear driftOSI driftDSI
for cell_n=1:n
    %for cell_n = [24]
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
        
        d = drift_all(cell_n);
        dr_spont(cell_n)=d.spont(1);
        tfs = squeeze(nanmean(d.orient_tune(1:2,:,:),3));
        tf_ratio(cell_n,:) = (tfs(2,:)-tfs(1,:))./(tfs(2,:) + tfs(1,:));
        %figure('Name',sprintf('cell %d',cell_n))
        color = {'b','r'};
        [sf_amp(cell_n) peak_sf(cell_n)] = max(squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1)));
        drift_amp(cell_n) = max(max(max(d.orient_tune)));
        for r=1:2  %%% tempfreq
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
        driftOSI(cell_n,:) = dr_osi(cell_n,pref_tf,:)
        driftDSI(cell_n,:) = dr_dsi(cell_n,pref_tf,:)
        driftF1F0(cell_n)=sum(d.orient_tune(pref_tf,2,:)) / sum(d.orient_tune(pref_tf,1,:));
    end
end

thresh=3;
mv_osi(mv_amp<thresh)=nan;
mv_dsi(mv_amp<thresh)=nan;
mv_spd(mv_amp<thresh)=nan;
mv_spd_ratio(mv_amp<thresh)=nan;

n_obs=n_obs+1; obs_name{n_obs} = 'mv osi'; obs(:,n_obs) =mv_osi;
n_obs=n_obs+1; obs_name{n_obs} = 'mv dsi'; obs(:,n_obs) = mv_dsi;
n_obs=n_obs+1; obs_name{n_obs} = 'mv spd'; obs(:,n_obs) = mv_spd;
n_obs=n_obs+1; obs_name{n_obs} = 'mv sp ratio'; obs(:,n_obs) = mv_spd_ratio;

d_thresh=1;
tf_ratio(tf_ratio>1)=1;
tf_ratio(tf_ratio<-1)=-1;
tf_ratio(sf_amp<d_thresh)=nan;
peak_sf(sf_amp<d_thresh)=nan;
driftDSI(sf_amp<d_thresh,:)=nan;
driftOSI(sf_amp<d_thresh,:)=nan;
driftF1F0(sf_amp<d_thresh)=nan;

[labels cmap clim]= makeColors(driftOSI(:,1),nan,'seq','GnBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels,5:11);
colormap(cmap); colorbar; set(gca,'Clim',clim); title('drift OSI');

[labels cmap clim]= makeColors(driftDSI(:,1),nan,'seq','GnBu');
f = plotSections(sections,anatomy,histox,histoy,histSection,labels,5:11);
colormap(cmap); colorbar; set(gca,'Clim',clim); title('drift DSI');

n_obs=n_obs+1; obs_name{n_obs} = 'tf0ratio'; obs(:,n_obs) = tf_ratio(:,1)';
n_obs=n_obs+1; obs_name{n_obs} = 'peakSF'; obs(:,n_obs) = peak_sf;
%n_obs=n_obs+1; obs_name{n_obs} = 'drift amp'; obs(:,n_obs) = sf_amp;

n_obs=n_obs+1; obs_name{n_obs} = 'driftDSI0'; obs(:,n_obs) = driftDSI(:,1)';
n_obs=n_obs+1; obs_name{n_obs} = 'driftDSI1'; obs(:,n_obs) = driftDSI(:,2)';
n_obs=n_obs+1; obs_name{n_obs} = 'driftOSI0'; obs(:,n_obs) = driftOSI(:,1)';
n_obs=n_obs+1; obs_name{n_obs} = 'driftOSI1'; obs(:,n_obs) = driftOSI(:,2)';
n_obs=n_obs+1; obs_name{n_obs} = 'dr spont'; obs(:,n_obs) = dr_spont;
n_obs=n_obs+1; obs_name{n_obs} = 'dr f1f0'; obs(:,n_obs) = driftF1F0;

n_obs=n_obs+1; obs_name{n_obs} = 'wn_cr'; obs(:,n_obs) = wn_cr_dom;
n_obs=n_obs+1; obs_name{n_obs} = 'wn_cr_norm'; obs(:,n_obs)=wn_cr_norm;
n_obs=n_obs+1; obs_name{n_obs} = 'sbc'; obs(:,n_obs)=(wn_cr_dom<-2*10^3);
n_obs=n_obs+1; obs_name{n_obs} = 'wn adapt'; obs(:,n_obs)=wn_adapt;


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
hist(fl_sust(fl_amp>4 & fl_lag==0));


figure
hist(tmin./tmax)

figure
plot(fl_lat(fl_amp>2 & fl_lag==0),fl_sust(fl_amp>2 & fl_lag==0),'o' );


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

% for cell_n = find(fl_resp)
%     f= fl_all(cell_n);
%     if ~isempty(f.N)
%     figure
%     for r = 1:2
%         hold on
%         if r==1
%             plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)))
%         else
%             plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'r')
%         end
%     end
%     end
% end




figure
for i = 1:n_obs;
    subplot(5,7,i)
    [h b]=hist(obs(:,i));
    bar(b,h);
end

mn = mean(wvform(find(inside),:));
[wv_pc score wv_eig] = princomp(wvform(find(inside),:)-repmat(mn,sum(inside),1));
figure
plot(wv_eig);
figure
plot(wv_pc(:,1:5))
figure
plot(score(:,1:5))
figure
plot(score(:,1),score(:,2),'o')

%usedcells= find(inside & ~(isnan(obs(:,4)) & isnan(obs(:,16))&isnan(obs(:,24))));
resp =~isnan(obs(:,find(strcmp(obs_name,'sta amp'))))+  ~isnan(obs(:,find(strcmp(obs_name,'fl_type')))) +...
    ~isnan(obs(:,find(strcmp(obs_name,'mv spd'))))+  ~isnan(obs(:,find(strcmp(obs_name,'dr f1f0')))) ;
    
sbc = wn_cr_dom<-2*10^3;
sbc=sbc';

usedcells = find(inside & resp>=1| sbc)

obs_used = obs(usedcells,:);
obs_used=replaceNan(obs_used,obs_name,'driftDSI0',0);
obs_used=replaceNan(obs_used,obs_name,'driftOSI0',0);
obs_used=replaceNan(obs_used,obs_name,'mv osi',0);
obs_used=replaceNan(obs_used,obs_name,'mv dsi',0);
obs_used=replaceNan(obs_used,obs_name,'sta amp',0);
obs_used=replaceNan(obs_used,obs_name,'mv spd',3);
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


c = nancov(obs_norm);
figure
imagesc(c)

cc = diag(c);
denom = sqrt(repmat(cc,1,length(cc)).*repmat(cc',length(cc),1));
nancorr = c./denom;

figure
imagesc(nancorr);
set(gca,'YTickLabel',labels)
set(gca,'YTick',1:length(obs_name))
title('nancor')

labels  ={};
for i = 1:length(obs_name);
    labels{i} = sprintf('%s  %d',obs_name{i},i);
end

figure
pcorr = corrcoef(obs_norm,'rows','pairwise');
imagesc(pcorr)
set(gca,'YTickLabel',labels)
set(gca,'YTick',1:length(obs_name))
title('pairwise')

labels  ={};
for i = 1:length(obs_name);
    labels{i} = sprintf('%s  %d',obs_name{i},i);
end

start_obs=7;
used_obs=start_obs:n_obs;
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'fl_lag'))); %%% fl_lag (1,2 are the only meaningful)
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wn_cr_norm'))); %%% cr_norm (captured by cr)
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'driftDSI1'))); %%% drift DSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'driftOSI1')));  %%% drift OSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wx'))); %%% drift DSI F1
used_obs=used_obs(used_obs~=find(strcmp(obs_name,'wy')));  %%% drift OSI F1


e = eig(pcorr(used_obs,used_obs));
figure
plot(e(end:-1:1));
[v d] = eig(pcorr(used_obs,used_obs));
figure
imagesc(v(:,size(v,2):-1:1));
set(gca,'YTickLabel',labels(used_obs))
set(gca,'YTick',1:length(obs_name))
title('pairwise corr')


v(abs(v)<0.1)=0;
% gamma = 1.5;
% v= sign(v).*abs(v).^gamma;
% figure
% imagesc(v(:,size(v,2):-1:1));
clear pc
for i=1:length(obs_norm);
    for j = 1:length(v);
        pc(i,j) = nansum(obs_norm(i,used_obs).*v(:,length(v)-j+1)');
        %pc(i,j) = nansum(obs_norm(i,6:end).*v(:,length(v)-j+1)')/sum(~isnan(obs_norm(i,6:end)));
        %pc(i,j) = nansum(obs_norm(i,6:end).*v(:,length(v)-j+1)')/sqrt(nansum(obs_norm(i,6:end).^2));
    end
end

figure
plot3(pc(:,1),pc(:,2),pc(:,3),'o')

% Y = pdist(pc(:,1:6), 'euclid');
%        Z = linkage(Y, 'ward');
%        T = cluster(Z, 'maxclust', n_clust);
% figure
% dendrogram(Z)
%%% n=7, pc1:7 phi=1.2

clear DB CH Dunn KL Han
close all
for rep=1:6;
    for c_rep=1:6
       npca=2*rep+4; 
        %npca =6;
        n_clust=5+c_rep;
       % n_clust=c_rep+5;
        
        phi=1.1                       % fuzzy exponent
        
        maxiter=200;
        toldif=0.0000001;
        
        ndata = size(pc, 1);         % number of data
        % initialise random
        Uinit= initmember(0.1,n_clust,ndata);
        tic
        data=pc(:,1:npca);
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
        imagesc(obs_sort');
        set(gca,'YTickLabel',labels(1:n_obs))
        set(gca,'YTick',1:length(obs_name))
        figure
        imagesc(obs_sort_mean');
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

%for i = 1:max(T);
for i = 1:max(T)
    id = nan(size(types));
    id(~isnan(types))=0;
    id(types==i)=1;
    
    [labels cmap clim]= makeColors(id,nan,'qual','Set1');
    f = plotSections(sections,anatomy,histox,histoy,histSection,labels);
    colormap(cmap); colorbar; set(gca,'Clim',clim); title(sprintf('type %d',i));
end



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
        title(sprintf('tp=%d sust=%0.1f c=%0.1f',fl_type(cell_n),fl_sust(cell_n),fl_onoffcorr(cell_n)));
        
        subplot(5,4,11);
        plot(fl_sztune(i,:));
        xlim([1 6])
        set(gca,'XTick',[]); xlabel('size');
    end
    
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
    
    
    subplot(5,4,4)
    plot(sections(histSection(cell_n)).coords(:,1),sections(histSection(cell_n)).coords(:,2),'b.'); hold on
    plot(histox(cell_n),histoy(cell_n),'r*');
    hold on
    plot(alignedX(cell_n),alignedY(cell_n),'g*')
    title(sprintf('Pax AP = %d',histSection(cell_n)));
    axis([-500 500 -500 500]);
    axis off
    
    subplot(5,4,8);
    lgnxy = anatomy(site(cell_n)).LGN;
    plot(lgnxy(:,2),lgnxy(:,1),'b.'); hold on
    plot(histox(cell_n),histoy(cell_n),'r*');
    axis([-500 500 -500 500]);
    axis off
    title(sprintf('monitor %0.0f',offsetX(cell_n)));
    
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
            title(sprintf('tf %0.1f os %0.1f ds %0.1f',tf_ratio(cell_n,f),driftOSI(cell_n,f),driftDSI(cell_n,f)));
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





