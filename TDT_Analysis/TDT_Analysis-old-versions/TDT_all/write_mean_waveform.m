%function create_analysisfile
%%%writes out the waveform/cluster data for selected units
%%%into an excel worksheet

clear all
close all
pack
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;

%%% cells = [1 5; 5 3;  9 2; 9 4; 9 8;  13 4; 13 5; 13 11]  % for 060607

cells = [1 4; 1 7; 1 8; 1 10; 5 2;  9 1; 9 4; 9 6; 13 1; 13 2; 13 5; 13 7; 13 8];

cells = [ 9 5; 9 9; 9 11; 13 2; 13 7]
[afname apname] = uigetfile('*.mat','analysis file');

load(fullfile(apname, afname))

clusterfilename = fullfile(pname,fname);
clear nspikes L_ratio trough_depth peak_height trough_width trough2peak

% %% merge two clusters
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
     c = cells(cell,2);

    
    figure
    plot(mean_wvform(:,cells(cell,1):cells(cell,1)+3,c));
     title(sprintf('ch %d cl %d',ch,c));
     
     [y mainchan]= min(mean_wvform(6,cells(cell,1):cells(cell,1)+3,c))
     wv(:,cell) = mean_wvform(:,mainchan+cells(cell,1)-1,c);
     wv(:,cell) = wv(:,cell)/abs(min(wv(:,cell)));
     

end
xlswrite('fullwvform',wv'); 
figure
plot(wv)
