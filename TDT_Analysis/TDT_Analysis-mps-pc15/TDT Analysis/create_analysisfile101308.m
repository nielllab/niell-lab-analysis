%function create_analysisfile
%%%writes out the waveform/cluster data for selected units
%%%into an excel worksheet

clear all
pack
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;

cells = [1 4 ; 1 7; 1 8; 1 9; 5 7; 9 3; 13 2; 13 4; 13 9; 13 10; 13 12]

[afname apname] = uiputfile('*.mat','analysis file');

fullaname = fullfile(apname, afname)

clusterfilename = fullfile(pname,fname);
clear nspikes L_ratio trough_depth peak_height trough_width trough2peak1

% %% merge two clusters = uncomment this as needed
% ch = 3;
% c1=1;
% c2=11;
% idx=idx_all(4*(ch-1)+1,:);
% idx(idx==c2)=c1;
% idx_all(4*(ch-1)+1,:)=idx;
% csize(ch,c1) = csize(ch,c1)+csize(ch,c2);
% Lratio(ch,c1) = 0.5*(Lratio(ch,c1)+Lratio(ch,c2));
% trough(ch,c1) = 0.5*(trough(ch,c1)+trough(ch,c2));
% peak(ch,c1) = 0.5*(peak(ch,c1)+peak(ch,c2));
% half_width(ch,c1) = 0.5*(half_width(ch,c1)+half_width(ch,c2));
% peak_trough(ch,c1) = 0.5*(peak_trough(ch,c1)+peak_trough(ch,c2));
% mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) = 0.5*(mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) + ...
%     mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c2));
% %% done merge

tet_data=1;
clear wv

for cell = 1:size(cells,1);
    if tet_data
        ch = ceil(cells(cell,1)/4);
    else ch = cells(cell,1);
    end
    ch
    
    %%% waveform and cluster parameters of interest
    c = cells(cell,2);
    nspikes(cell) = csize(ch,c);
    L_ratio(cell) = Lratio(ch,c);
    trough_depth(cell) = trough(ch,c);
    peak_height(cell) = peak(ch,c);
    trough_width(cell) = half_width(ch,c);
    trough2peak(cell) = peak_trough(ch,c);

    figure
    plot(mean_wvform(:,cells(cell,1):cells(cell,1)+3,c));
    title(sprintf('ch %d cl %d',ch,c));

    [y mainchan]= min(mean_wvform(6,cells(cell,1):cells(cell,1)+3,c))
    wv(:,cell) = mean_wvform(:,mainchan+cells(cell,1)-1,c);
    wv(:,cell) = wv(:,cell)/abs(min(wv(:,cell)));
 
      %%% calculate ISI
    ch = cells(cell,1);
    t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == c)));  
    dt = diff(t1);
    dt =dt(dt<.02);
    figure
    hist(dt,.0005:0.001:.02);
    
    %%% calculate autocorrelation 
    %%% to do this fast, do ISI for neighboring spikes, then spikes
    %%% separated by 1, 2, 3 ... up to 12
    for s =1:12;
        if s==1
            dt = diff(t1,s);
            dt= dt(dt<.03);
            h = hist(dt,.0005:0.001:.03);
        else
            dt= t1(s:length(t1))- t1(1:length(t1)-s+1);
            dt=dt(dt<.03);
            h=h+hist(dt,.0005:0.001:.03);
        end
    end
    figure
    bar(h);
    title(sprintf('ch %d cl %d',ch,c));
end
xlswrite('fullwvform',wv');
figure
plot(wv)

save(fullaname, 'clusterfilename', 'Tank_Name',  'cells', 'nspikes', 'L_ratio', 'trough_depth', 'peak_height', 'trough_width', 'trough2peak','idx_all');

