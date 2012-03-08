function create_analysisfile
%%%writes out the waveform/cluster data for selected units
%%%into an excel worksheet

global goodcells;

[fname, pname] = uigetfile('*.mat','cluster data')
oldpname = pname;
oldfname = fname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;

[basename] =sscanf(oldfname,'cluster_data_%s.mat')
basename = basename(1:length(basename)-4)

%%% cells = [1 5; 5 3;  9 2; 9 4; 9 8;  13 4; 13 5; 13 11]  % for 060607

[afname apname] = uiputfile('*.mat','analysis file');

fullaname = fullfile(apname, afname)

clusterfilename = fullfile(pname,fname);
clear nspikes L_ratio trough_depth peak_height trough_width trough2peak1

tet_data=1;
clear wv

linecolor = [0 0 1; 0 1 0 ; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; .25 0 0.5; 0 0.5 0 ; 0.5 .25 0; 0.5 0 1; 0 0.5 0.5];
ncells = 0;
for tet=1:4
    close all
    tet_ch = (tet-1)*4 +1;
    goodcells = zeros(1,12);
    try
        open(sprintf('%shist%s_t%d.fig',oldpname,basename,tet))
        set(gcf,'Position',[10 50 500 400]);
        open(sprintf('%ssnip%s_t%d.fig',oldpname,basename,tet))
        set(gcf,'Position',[10 550 500 400]);
        open(sprintf('%sclust%s_t%d.fig',oldpname,basename,tet))
        set(gcf,'Position',[600 50 500 400]);
    end
    select_fig = figure;
    set(gcf,'Position',[600 600 500 400]);
    for i =1:12;
        subplot(3,4,i);

        plot([0 1],[0 0],'LineWidth',6,'Color',linecolor(i,:));
        axis off
        op = get(gca,'OuterPosition');
        axes('position',op); axis off;
        a = patch([0 0 1 1],[0 1 1 0],'w');
        set(a,'FaceAlpha',0)
        set(a,'EdgeAlpha',0)
        set(a,'ButtonDownFcn',@togglegoodcell);
        set(a,'UserData',i); % store cell number in axes
    end
    done=0;
    while ~done
        k = waitforbuttonpress;
        if k==1 %%% keyboard
            done=1;
        end
    end
    used = find(goodcells);

    for i =length(used):-1:1;

        c = used(i);

        figure
        subplot(1,2,1);
        plot(mean_wvform(:,tet_ch:tet_ch+3,c));
        title(sprintf('ch %d cl %d',tet,c));
      

       subplot(1,2,2);
        t1= squeeze(event_times_all(tet_ch,find(idx_all(tet_ch,:) == c)));
        dt = diff(t1);
        dt =dt(dt<.02);
        hist(dt,.0005:0.001:.02);
          set(gcf,'Position',[50 50 800 400], 'Color',linecolor(c,:));
    end

    figure(select_fig);
    done=0;
    while ~done
        k = waitforbuttonpress;
        if k==1 %%% keyboard
            done=1;
        end
    end

    used = find(goodcells);

    cells(ncells+1:ncells+length(used),1)=tet_ch;
    cells(ncells+1:ncells+length(used),2)=used;
    ncells = ncells+length(used);
end

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
    wv(:,cell) = mean_wvform(:,mainchan+cells(cell,1)-1,c);
    wv(:,cell) = wv(:,cell)/abs(min(wv(:,cell)));

end
xlswrite('fullwvform',wv');
figure
plot(wv)


save(fullaname, 'clusterfilename', 'Tank_Name',  'cells', 'nspikes', 'L_ratio', 'trough_depth', 'peak_height', 'trough_width', 'trough2peak','idx_all');

