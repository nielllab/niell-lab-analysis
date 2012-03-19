clear all
afile = {'C:\data\ephys matlab data\030912_LGN\recording3\analysis3.mat' ...
    'C:\data\ephys matlab data\031512_LGN\analysis2.mat' ...
    'C:\data\ephys matlab data\031612_lgn\recording1\analysis1.mat'}
n=0;
for i = 1:length(afile);
    i
    load(afile{i});
    cell_range = n+(1:length(drift));
    n=n+length(drift);
    if size(wn,2)==2
        wn_all(cell_range,:)=wn;
    else
        wn_all(cell_range,1)=wn;
    end
    drift_all(cell_range)=drift;
    mv_all(cell_range)=mv;
    fl_all(cell_range)=fl;
    site(cell_range)=i;
    cell_id(cell_range,:) = cells;
end

clear p
close all
n_sta=0;

for cell_n = 1:n
    cell_n
    wn_all(cell_n).sta_fit=[nan nan nan nan nan nan ];
    for s =1:3
        sta = double(squeeze(wn_all(cell_n).svd_xy(s,:,:)));
        
        [fit g]= fitLGNrf(sta);
        background = find(abs(g-fit(2))<(0.1*abs(fit(1))));
        if fit(1)>1 | fit(3)<1 | fit(4)<1 | fit(3)>size(sta,1) | fit(4)>size(sta,2)
            break
        end
        z=abs(fit(1))/std(sta(background))
        
        if z>6
            wn_all(cell_n).sta_t = wn_all(cell_n).svd_t(:,s);
            if wn_all(cell_n).sta_t(6)<0;
                wn_all(cell_n).sta_t = wn_all(cell_n).sta_t *-1;
                fit(1) = fit(1)*-1;
                fit(2)=fit(2)*-1;
            end
            wn_all(cell_n).sta_fit=fit;
            wn_all(cell_n).sta_final=sta;
            wn_all(cell_n).zscore=z;
%             figure
%             subplot(2,2,1)
%             imagesc(sta); axis equal; axis tight;
%             
%             subplot(2,2,2)
%             imagesc(g);axis equal; axis tight;
%             
%             subplot(2,2,4);
%             plot(wn_all(cell_n).sta_t)
%             xlim([0 30])
%             
            
            n_sta=n_sta+1;
            break
        end
    end
end

ind=0;
tmin=zeros(1,cell_n); tmax = tmin; wn_lat=tmin;

for cell_n=1:n
    wn_resp(cell_n) = wn_all(cell_n).responsiveness(1);
    wn_phase(cell_n) =wn_all(cell_n).phase;
    A(cell_n) = wn_all(cell_n,1).sta_fit(1);
    wx(cell_n)=wn_all(cell_n,1).sta_fit(6);
    wy(cell_n)=wn_all(cell_n,1).sta_fit(5);
    x0(cell_n)= wn_all(cell_n,1).sta_fit(4);
    y0(cell_n)= wn_all(cell_n,1).sta_fit(3);
    if ~ isempty(wn_all(cell_n,1).sta_t)
        [z wn_lat(cell_n)] = max(wn_all(cell_n,1).sta_t);
        tmax(cell_n)=max(wn_all(cell_n,1).sta_t)*sign(A(cell_n));
        tmin(cell_n)=min(wn_all(cell_n,1).sta_t)*sign(A(cell_n));
    end
    wn_cr(cell_n) = (wn_all(cell_n).crf(10) - wn_all(cell_n).crf(1));
    wn_cr_norm(cell_n) = (wn_all(cell_n).crf(10) - wn_all(cell_n).crf(1)) ./ ( (wn_all(cell_n).crf(10) + wn_all(cell_n).crf(1)));
end

for cell_n = 1:n
    f= fl_all(cell_n);
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
    
    use_sizes = find(sz>0.5*max(sz))
    for r=1:2
        
        h = squeeze(mean(f.onset_hist(r,use_sizes,:),2))-f.spont;
        fl_tseries(cell_n,r,:) = h;
        [onset(cell_n,r) onset_lat(cell_n,r)]=max(h(f.onset_bins>=0.3 & f.onset_bins<0.55));
        onset_trans(cell_n,r) = mean(h(f.onset_bins>=0.3 & f.onset_bins<0.55))/onset(cell_n,r);
        [offset(cell_n,r) offset_lat(cell_n,r)]=max(h(f.onset_bins>=0.55 & f.onset_bins<0.8));
        offset_trans(cell_n,r) = mean(h(f.onset_bins>=0.55 & f.onset_bins<0.8))/offset(cell_n,r);
    end
    
    fl_lat(cell_n)=onset_lat(cell_n,fl_onoff(cell_n));
    fl_trans(cell_n) =onset_trans(cell_n,fl_onoff(cell_n));
    
    fl_lag(cell_n) = f.lag;
    fl_onoff(cell_n) = 2*(fl_onoff(cell_n)-1.5);
    if fl_lag(cell_n)==1  %%% if a cell responds with 1 frame delay, it is responding to offset
        fl_onoff(cell_n)=-1*fl_onoff(cell_n);
    end
end


for cell_n=1:n
%for cell_n = [24]
    m= mv_all(cell_n);
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
    use_sz = find(sz>0.5*max(sz))
    sp = squeeze(mean(nanmean(filled_h,3),1))-m.spont;
    use_sp = find(sp>0.5*max(sp))
    os_tune = squeeze(mean(nanmean(filled_h(use_sz,use_sp,:),2),1))-m.spont
    os_N= squeeze(sum(nansum(m.n_unique(onoff,use_sz,use_sp,:),3),2))
  
    [mv_amp(cell_n) mv_sz(cell_n)] =max(sz);
    if length(sz)==5
        mv_sz(cell_n)=mv_sz(cell_n)-1;
    end
    [empt mv_spd(cell_n)] = max(sp);
    [mv_osi(cell_n) mv_prefO(cell_n)]= calcOSI(os_tune,0);
    [mv_dsi(cell_n) mv_prefD(cell_n)]= calcOSI(os_tune,1)
    
       figure('Name',sprintf('cell %d',cell_n))
    subplot(2,2,1)
    plot(sz);
    axis([1 3 min(0,min(sz)) max(1,max(sz))])
    subplot(2,2,2)
    plot(sp);
    axis([1 5 min(0,min(sp)) max(1,max(sp))])
    subplot(2,2,3)
    plot(os_tune)
    axis([1 8 min(0,min(os_tune)) max(1,max(os_tune))])
    title(sprintf('OSI = %f DSI=%f',mv_osi(cell_n),mv_dsi(cell_n)))
    
    d = drift_all(cell_n);
    tfs = squeeze(nanmean(d.orient_tune(1:2,:,:),3));
    tf_ratio(cell_n,:) = (tfs(2,:)-tfs(1,:))./(tfs(2,:) + tfs(1,:));
    figure('Name',sprintf('cell %d',cell_n))
    color = {'b','r'}
    [sf_amp peak_sf(cell_n)] = max(squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1)));
    
    for r=1:2
        for f=1:2
            [peak_amp peak_tf(cell_n)] = max(squeeze(nanmean(d.sf_tune(:,f,:),3)));
            drift_os = squeeze(d.orient_tune(r,f,:));
            
            subplot(2,2,f)
            hold on
            plot(drift_os,color{r});
            title(sprintf('tfratio %f',tf_ratio(cell_n,f)));
            
            drift_sf= squeeze(d.sf_tune(r,f,:));
            subplot(2,2,f+2)
            hold on
            plot(drift_sf,color{r});
            title(sprintf('peak sf %d',peak_sf(cell_n)));
            %axis([1 8 min(0,min(drift_os)) max(drift_os)])
        end
    end
end
    
    
    
    
    
    

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

%%% wn on/off position
figure
for i =1:2
    hold on
    if i==1
        plot(x0(A<0),y0(A<0),'o');
    elseif i==2
        plot(x0(A>0),y0(A>0),'go');
    end
end
axis ij


figure
plot(tmin./tmax,wx,'o')

figure
plot(fl_lag,wx,'o')
figure
plot(x0,wx,'o')
axis([0 128 0 15])
figure
plot(x0,wy,'o')

figure
hist(fl_lat(fl_amp>4 & fl_lag==0));

figure
hist(fl_trans(fl_amp>4 & fl_lag==0));


figure
hist(tmin./tmax)

figure
plot(fl_lat(fl_amp>2 & fl_lag==0),fl_trans(fl_amp>2 & fl_lag==0),'o' );


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
axis ij

fl_responder=find(onset>5)
figure
plot(onset_trans(fl_responder),onset_lat(fl_responder),'o')

figure
hist(onset_trans(onset>5));

figure
hist(onset_lat(onset>5),0:10)

figure
plot(fl_sz(fl_amp>1),fl_amp(fl_amp>1),'o')

figure
hist(fl_supp(fl_amp>2))

figure
plot(fl_lag(fl_amp>4),fl_supp(fl_amp>4),'o')


fl_resp = find(fl_amp>4 | (fl_amp>2 & fl_spont<1))

for cell_n = find(fl_resp)
    f= fl_all(cell_n);
    
    figure
    for r = 1:2
        hold on
        if r==1
            plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)))
        else
            plot(f.onset_bins(1:end-1),squeeze(mean(f.onset_hist(r,:,:),2)),'r')
        end
    end
end