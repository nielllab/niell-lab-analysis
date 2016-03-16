function select_units
%%% after clustering is done, this program lets you go through and select
%%% the clusters to use for analysis
%%% written by Cris Niell, 2009-2010
psfilename = 'C:\tempPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end

dbstop if error
global goodcells;

%%% read in the results of clustering
[fname, pname] = uigetfile('*.mat','cluster data')
oldpname = pname;
oldfname = fname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;


basename = oldfname(14:end-4)

% [basename] =sscanf(oldfname,'cluster_data_%s.mat')
% basename = basename(1:length(basename)-4)

%%% where to save analysis results
[afname apname] = uiputfile('*.mat','analysis file');

fullaname = fullfile(apname, afname)

clusterfilename = fullfile(pname,fname);
clear nspikes L_ratio trough_depth peak_height trough_width trough2peak1

tet_data=1;
clear wv

% %%% merge two clusters
% %%% this is a total hack! but can't think of a good gui to do this
% ch = 3;   %%% ch is tetrode number, i.. 1-16
% c1=2;
% c2=4;
% idx=idx_all{4*(ch-1)+1};
% idx(idx==c2)=c1;
% idx_all{4*(ch-1)+1}=idx;
% csize(ch,c1) = csize(ch,c1)+csize(ch,c2);
% Lratio(ch,c1) = 0.5*(Lratio(ch,c1)+Lratio(ch,c2));
% trough(ch,c1) = 0.5*(trough(ch,c1)+trough(ch,c2));
% peak(ch,c1) = 0.5*(peak(ch,c1)+peak(ch,c2));
% half_width(ch,c1) = 0.5*(half_width(ch,c1)+half_width(ch,c2));
% peak_trough(ch,c1) = 0.5*(peak_trough(ch,c1)+peak_trough(ch,c2));
% mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) = 0.5*(mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) + ...
%     mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c2)); 
% % %%% done merge
% 
% %%% merge two clusters
% %%% this is a total hack! but can't think of a good gui to do this
% ch = 4;   %%% ch is tetrode number, i.. 1-16
% c1=2;
% c2=3;
% idx=idx_all{4*(ch-1)+1};
% idx(idx==c2)=c1;
% idx_all{4*(ch-1)+1}=idx;
% csize(ch,c1) = csize(ch,c1)+csize(ch,c2);
% Lratio(ch,c1) = 0.5*(Lratio(ch,c1)+Lratio(ch,c2));
% trough(ch,c1) = 0.5*(trough(ch,c1)+trough(ch,c2));
% peak(ch,c1) = 0.5*(peak(ch,c1)+peak(ch,c2));
% half_width(ch,c1) = 0.5*(half_width(ch,c1)+half_width(ch,c2));
% peak_trough(ch,c1) = 0.5*(peak_trough(ch,c1)+peak_trough(ch,c2));
% mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) = 0.5*(mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c1) + ...
%     mean_wvform(:,4*(ch-1)+1:4*(ch-1)+4,c2)); 
% % %%% done merge

linecolor = [0 0 1; 0 1 0 ; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; .25 0 0.5; 0 0.5 0 ; 0.5 .25 0; 0.5 0 1; 0 0.5 0.5];
ncells = 0;
ceil(length(idx_all)/4)

for tet=1:ceil(length(idx_all)/4);   %%% for each tetrode, show histograms, waveforms, and ica scatter plots
    close all
    tet_ch = (tet-1)*4 +1;
    goodcells = zeros(1,12);
    try
        open(sprintf('%shist%s_t%d.fig',oldpname,basename,tet))
        sprintf('%shist%s_t%d.fig',oldpname,basename,tet)
        set(gcf,'Position',[10 50 500 400]);
        open(sprintf('%ssnip%s_t%d.fig',oldpname,basename,tet))
        set(gcf,'Position',[10 550 500 400]);
        for i = 1:4
            subplot(2,2,i); ylim([-100 50])
        end
        
        open(sprintf('%sclust%s_t%d.fig',oldpname,basename,tet))
        set(gcf,'Position',[600 50 500 400]);
    end

    %%% make a toggle figure to allow selection of units
    select_fig = figure;
    set(gcf,'Position',[600 600 500 400]);
    for i =1:12;
        subplot(3,4,i);

        plot([0 1],[0 0],'LineWidth',6,'Color',linecolor(i,:));
        axis off
        op = get(gca,'OuterPosition');
        axes('position',op); axis off;
        a = patch([0 0 1 1],[0 1 1 0],'w','PickableParts','All');
        set(a,'FaceAlpha',0)
        set(a,'EdgeAlpha',0)
        set(a,'ButtonDownFcn',@togglegoodcell);
        set(a,'UserData',i); % store cell number in axes
        if i<=max(idx_all{tet_ch})
            title(sprintf('%0.2f', Lratio(tet,i)));
        end
    end
    
    pause
    used = find(goodcells);

    %%% once potential clusters have been selected, show avg waveform and ISI
   
    
    for i =length(used):-1:1;

        c = used(i);

        figure
        subplot(2,2,1);
        plot(mean_wvform(:,tet_ch:tet_ch+3,c));
        title(sprintf('ch %d cl %d',tet,c));
      
       subplot(2,2,2);
        t1= squeeze(event_times_all{tet_ch}(find(idx_all{tet_ch} == c)));
        length(t1)
        dt = diff(t1);
        dt =dt(dt<.02);
        hist(dt,.0005:0.001:.02);
          set(gcf,'Position',[50 50 800 400], 'Color',linecolor(c,:));
         
         if exist('wave_all','var')
             wvall = wave_all{tet};
           wvclust = wvall(find(idx_all{tet_ch}==c),:,:);
            
          dt = diff(t1);
          breaks = find(dt>5*10^4);
          breaks(end+1)=length(dt);
          dt(dt>5*10^4)=1;
          tmerge = cumsum(dt);
                binwidth=60;
                amps =squeeze(min(wvclust(:,5:10,:),[],2));
                 clear ampmean
                 for t =1:floor(max(tmerge)/binwidth);
                      ampmean(t,:) = median(amps(tmerge>(t-1)*binwidth & tmerge<t*binwidth,:),1);
                  end
                  
        subplot(2,2,3:4)
      if length(amps)~=0  %%% happens from merge
          plot([0 tmerge],amps,'.','MarkerSize',2 ); hold on
        plot((binwidth:binwidth:max(tmerge))-binwidth/2,ampmean,'LineWidth',2);
        for b = 1:length(breaks)
            plot([tmerge(breaks(b)) tmerge(breaks(b))],[-10^-4 2*10^-5])
            bl = Block_Name{b}; bl(bl=='_')=' ';
            text(tmerge(breaks(b))-120,2*10^-5,bl,'Rotation',90,'HorizontalAlignment','right')
        end
      end
      
        
         end
         
          %%% call Erik's code to calculate cluster separation
          %%%showClusterSeparation(wave_all{tet},idx_all{tet_ch},i);
          set(gcf, 'PaperPositionMode', 'auto');
          print('-dpsc',psfilename,'-append');
    end

    %%% give the user a chance to revise their choices
    figure(select_fig);
    pause

    used = find(goodcells);

    cells(ncells+1:ncells+length(used),1)=tet_ch;
    cells(ncells+1:ncells+length(used),2)=used;
    ncells = ncells+length(used);
    close all
end

clear wv
cells
tet_data=1;
%%% consolidate data for the toggle cells selected
for cell = 1:size(cells,1);
    if tet_data
        ch = ceil(cells(cell,1)/4);
    else ch = cells(cell,1);
    end
    ch
    c = cells(cell,2);
    nspikes(cell) = csize(ch,c);
    L_ratio(cell) = Lratio(ch,c);
    trough_depth(cell) = trough(ch,c);
    %trough_depth(cell) = NaN;
    peak_height(cell) = peak(ch,c);
    trough_width(cell) = half_width(ch,c);
    trough2peak(cell) = peak_trough(ch,c);

    [y mainchan]= min(mean_wvform(6,cells(cell,1):cells(cell,1)+3,c))
    peakchan(:,cell) = (ch-1)*4 + mainchan;
    wv(:,cell) = mean_wvform(:,mainchan+cells(cell,1)-1,c);
    wv(:,cell) = wv(:,cell)/abs(min(wv(:,cell)));
    spikeT{cell} = event_times_all{cells(cell,1)}(find(idx_all{cells(cell,1)} == c));
end

figure
plot(wv)
cells
peakchan

save(fullaname, 'clusterfilename', 'Tank_Name',  'cells', 'nspikes', 'L_ratio', 'trough_depth', 'peak_height', 'trough_width', 'trough2peak','idx_all', 'pname','wv','spikeT','peakchan');

[f p] = uiputfile('*.pdf','save pdf');
if f~=0
    try
   ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
catch
    display('couldnt generate pdf');
    end
end
delete(psfilename);
