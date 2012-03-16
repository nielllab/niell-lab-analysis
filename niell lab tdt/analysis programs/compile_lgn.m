clear all
afile = {'C:\data\ephys matlab data\030912_LGN\recording3\analysis3.mat' ...
    'C:\data\ephys matlab data\031512_LGN\analysis2.mat'}
n=0;
for i = 1:length(afile);
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


                n_sta=n_sta+1;
            break
        end
    end
end

ind=0;
tmin=zeros(1,cell_n); tmax = tmin; wn_lat=tmin;

for cell_n=1:n
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
        
end

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
axis equal

%%% wn latency
figure
hist(wn_lat(wn_lat>0))

%%% wn sites
figure
color = {'bo','go','ro','co','mo','yo','ko'; 'bx','gx','rx','cx','mx','yx','kx'};
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

for cell_n = 1:n
   f= fl_all(cell_n);
%     flash_spont = mean(f.onset_hist(f.onset_bins<0.25))
%     figure
%     plot(f.onset_hist')
   [junk fl_onoff(cell_n) ] = max(mean(f.hist_all,2));
  
    sz= squeeze(f.hist_all(fl_onoff(cell_n),:))-f.spont;
    [fl_amp(cell_n) fl_sz(cell_n)] = max(sz) ;
    fl_supp(cell_n) = sz(end)/max(sz);
    fl_spont(cell_n) = f.spont;
   
    fl_lag(cell_n) = f.lag;
%      fl_onoff(cell_n) = 2*(fl_onoff(cell_n)-1.5);
%      if fl_lag(cell_n)==1  %%% if a cell responds with 1 frame delay, it is responding to offset
%          fl_onoff(cell_n)=fl_onoff(cell_n);
%      end
end


figure
plot(fl_amp(fl_amp>1),fl_sz(fl_amp>1),'o')
figure
hist(fl_supp(fl_amp>2))

for i = 1:2
figure
plot(fl_onoff(fl_lag==0),A(fl_lag==0),'o')
hold on
plot(fl_onoff(fl_lag==1),A(fl_lag==1),'ro')

end