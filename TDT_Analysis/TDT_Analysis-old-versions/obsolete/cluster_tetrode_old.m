%function cluster_tetrode
% Matlab codes for clustering tetrode data
% Reads in waveform snippets from TDT tank, performs alignment,
% calculates independent components, calls Klustakwik, then displays
% cluster parameters
% User can then use select_units.m to choose appropriate clusters
% written by Cris Niell, 2006-2010

close all
Tank_Name=[];

global goodcells;
done=0
nblock=0;

% option to redo an old clustering
recluster= input('reclustering? (y/n) ','s');
if recluster~='y'
    while ~done
        if nblock==0
            pname = uigetdir('C:\data\TDT','block data')   %%% start location of tanks for search
        else
            pname = uigetdir(selected_path,'block data')
        end
        if pname==0;
            done=1;
        else
            nblock=nblock+1;
            delims = strfind(pname,'\');
            selected_path = pname(1 :delims(length(delims))-1)
            Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
            Block_Name{nblock} = pname(delims(length(delims))+1 :length(pname))
        end
    end

    nblocks = size(Block_Name,2);
    TTX = openTTX(Tank_Name,char(Block_Name(1))); % to initialize

    Event_Name_Snip='Snip'
    trgname = 'Xtrg';
    max_time = 60*60; %%% total recording duration (secs)

    use_tets = 1:4;  %%% use all tetrodes
    badsites = [];   %%% no bad sites to be removed
else    %% reclustering
    [fname pname] = uigetfile('','cluster file');
    load(fullfile(pname,fname));
    use_tets = input('enter tetrodes to recluster, e.g. [1 4] : ');
    badsites = input('enter channel to zero out, e.g. [3 8 14] [](none) : ');
end


for tet=use_tets
    
    %%%% read in all waveforms for this tetrode
    tet_title = sprintf('tetrode %d',tet);
    clear X
    for tet_ch=1:4
        N=0;
        Nblock=0;
        for block= 1:nblocks;
            TTX = openTTX(Tank_Name,char(Block_Name(block))); % to initialize
            ch = (tet-1)*4+tet_ch;
            invoke(TTX,'CreateEpocIndexing');
            MyEpocs = invoke(TTX, 'GetEpocsV', trgname, 0,0, 1000);
            Select_Duration(1) = 0;
            Nblock(block)= invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch, 0, Select_Duration(1),max_time,'ALL');
            W = invoke(TTX, 'ParseEvV', 0, Nblock(block));
            if ismember(ch, badsites)
                W=zeros(size(W));
            end
            event_times_all(ch,N+(1:Nblock(block))) = invoke(TTX, 'ParseEvInfoV', 0, Nblock(block), 6) +10^5 *(block-1);  %   3  = event codes
            noise(ch) = std(W(1,:))/(10^-6);
            sprintf('channel %d : noise = %f',ch,noise(ch))
            X(N+(1:Nblock(block)),1:size(W,1),tet_ch)=W';
            % blockID(ch,N+(1:Nblock(block)))=block;
            N=N+Nblock(block);
            NS(tet_ch)= N;
        end
    end
    clear W
    noise = squeeze(std(X(:,1,:),[],1));

    %%% waveform data is now stored in X
    %%% first dimension is waveform #, second dimension is timepoint
    %%% third dimension is channel # (1-4)
    
    %%% sometimes TDT loses spikes.
    %%% if one channel is dropped for a snippet, then all the following
    %%% snippets are offset.
    %%% so we need to find any events that only have 3channels worth of
    %%% data, and remove them. we do this by looking at event times
    %%% note: shouldn't be necessary for a more reliable acquisition system!    
    if std(NS)~=0 | sum(std(event_times_all((tet-1)*4+1:4,:)))>0
        sprintf('warning missing spikes')

        Xold = X;
        etimes_old = event_times_all;
        event_times_all((tet-1)*4+(1:4),:)=0;

        t234=etimes_old((tet-1)*4+(2:4),:);
        t234(:,max(NS))=Inf;
        lastmatch(1:3)=0;
        N=0;
        nf=0;
        tm= zeros(NS(1),3);
        for i=1:NS(1)

            t=etimes_old((tet-1)*4+1,i);
            failed=0;

            for j=1:3
                inc=1;
                %                 t234(j,lastmatch(j)+inc)
                %                 t
                while (t234(j,lastmatch(j)+inc)~=t)
                    if t234(j,lastmatch(j)+inc) > t
                        failed=1;
                        break
                    else inc = inc+1;
                    end
                end
                if failed==1

                    break;
                end
                tm(N+1,j)=lastmatch(j)+inc;
                lastmatch(j) = lastmatch(j)+inc;
                failed;
            end
            if failed==0
                N=N+1;
                X(N,:,1)=Xold(i,:,1);
                event_times_all((tet-1)*4+(1:4),N)=t;
                for j=1:3
                    X(N,:,j+1) = Xold(tm(N,j),:,j+1);
                end
            else
                nf=nf+1;
                failedN(nf)=i;
            end
        end
        X=X(1:N,:,:);
        clear newtimes

        newtimes(:,1)  = etimes_old((tet-1)*4+(2),tm(1:N,1));
        newtimes(:,2)  = etimes_old((tet-1)*4+(3),tm(1:N,2));
        newtimes(:,3)  = etimes_old((tet-1)*4+(4),tm(1:N,3));
        max(std(newtimes,[],2))
        i
        N
    end

    %%% cut out beginning of waveform before threshold, no information there
    X = X(:,6:30,:);

    clear Xshift Xold etimes_old newtimes;
    pack

    %%% begin alignment of spikewaveforms

    trigT = 6;  %%% timepoint to align on
    shiftrange=6;  %%% maximum shift allowed in alignement

    sz = 25-shiftrange;
    Xshift = zeros(N,sz,4);
   
    %%% remove snippets that have an initial voltage that is beyond the
    %%% usual noise range, as they are generally artifacts
  
    threshold_voltage = 1;  %%% use thresholding?
    v0thresh = 40*10^-6;    %%% threshold for initial point (a real spike will not have such a large deviation before beginning of spike
    vmax_thresh = 800*10^-6;  %%% threshold for largest excursion - only artifacts will have amplitude greater than this

    v0 = squeeze(max(abs(X(:,1:3,:)),[],3));
    v0 = squeeze(max(v0,[],2));
    vmax = squeeze(max(abs(X(:,:,:)),[],3));
    vmax = squeeze(max(vmax,[],2));
    
    if threshold_voltage
        used = find(v0< v0thresh & vmax <vmax_thresh);
    else
        used = find(v0>0);
    end
    figure
    hist(v0);

    %%% perform alignement based on minimum point, as averaged across channels
    for chan = 1:4
        [Y I] = min(squeeze(mean(X(:,trigT:trigT+shiftrange,:),3)),[],2);
        for i = 1:N
            Xshift(i,:,chan) = X(i,I(i) : (I(i)+sz-1),chan);
        end
    end

    figure
    hist(I);
    X= Xshift(used,:,:);

    %%% remove common mode across the four channels
    %%% this step is optional - ICA does a pretty good job of this anyways
    avg = mean(X,3);
    for chan = 1:4
        X(:,:,chan)=X(:,:,chan)-avg;
    end
    
    clear Xshift Y c_score;
    pack

    Xsize= size(X);

    clear Xshift allXdata chandata
    pack

    min_range = 5:10;
    max_range = 14:22;
    ica_fig=figure
    score = zeros(size(X,1),8);

    %%% calculate ICA components
    %%% note, some variable names reflect the fact that this used to be pca, not ica.
    
    choose_pca=1;   %%% select the componets to use, or use all
    
    if choose_pca
        n_pca=8;
    else
        n_pca =8;
    end
    mn = mean(X,3);
    pc_data = [X(:,:,1) X(:,:,2) X(:,:,3) X(:,:,4)];  %%% concatenate waveforms to send to ICA
    size(pc_data)
    [s coeff u] = fastica(pc_data'*10^6,'numOfIC',n_pca,'lastEig',n_pca,'stabilization','on','g','tanh','approach','symm');
    score = s';
    %%% components are in coeff, and the ICA value is in score
    
    
   %%% display ICA components, and choose ones to use
   for i = 1:n_pca
        subplot(ceil(n_pca/4),4,i);
        sz = size(X,2);
        plot(coeff(1:sz,i));
        hold on
        plot(coeff(sz+1:2*sz,i),'g');
        plot(coeff(2*sz+1:3*sz,i),'r');
        plot(coeff(3*sz+1:4*sz,i),'c');
        %plot(coeff(4*sz+1:5*sz,i),'k');

        % layer a (semi-)transparent patch over the axes to trap clicks
        op = get(gca,'OuterPosition');
        axes('position',op); axis off;
        a = patch([0 0 1 1],[0 1 1 0],'w');
        set(a,'FaceAlpha',0)
        set(a,'EdgeAlpha',0)
        set(a,'ButtonDownFcn',@togglegoodcell);
        set(a,'UserData',i); % store cell number in axes
    end

    clear pc_data  I mn
    pack
    clear pca_use
    figure(ica_fig)

    goodcells = zeros(1,n_pca);

    if choose_pca
        pause    %%% wait until the user presses a key, allowing time to click on components to use
       pca_use = find(goodcells);
        n_pca = length(pca_use);

        %% % in some versions of matlab, pause doesn't give control to the
        %% figure, so use this instead
%         done=0
%           while ~done
%             k = waitforbuttonpress;
%             if k==1 %%% keyboard
%                 done=1;
%             end
%         end


    else
        pca_use=1:n_pca;
    end
 
    %%% send selected ICA scores to klustakwik
    score = score(:,pca_use);
    fname='pca.fet.1'
    dlmwrite(fname,n_pca);
    dlmwrite(fname,score,'-append','delimiter',' ');
    sprintf('clustering data .....')
    tic; dos('C:\KlustaKwik_amd pca 1 -MinClusters 5 -MaxClusters 10 -MaxPossibleClusters 12 -Verbose 0 -nStarts 2','-echo'); toc
    
    %%% klustakwik returns a file with first value # of cluster, then cluster # for each snippet
    cl_name = 'pca.clu.1';
    idx = uint8(dlmread(cl_name));
    n_clust = idx(1)
    if sum(idx==1)<10   %%% first cluster is outliers, if there's only a few then just remove them
        n_clust=n_clust-1;
        idx = idx-1;
    end
    idx = idx(2:size(idx,1));
    size(idx_all)
    size(idx)
    size(used)
    
    for t = 1:4
        t
        idx_all(((tet-1)*4+t),used)=idx';
    end

    %%% add back the common mode that was subtracted off waveforms above
    for chan = 1:4
        X(:,:,chan)=X(:,:,chan)+avg;
    end

    save waveformdata X idx
    
    for c = 1:n_clust
        clust_size(c)=sum(idx==c);
    end
    clear group
    clust_size
    group = zeros(max(clust_size),size(X,2),4);


    %%% select waveforms and calculate mean for each cluster
    %%% simplifies a lot of the code for display below
    
    for c = 1:n_clust
        members = find(idx==c);
        group(1:size(members,1),1:size(X,2),1:4)=X(members,:,:);
        mean_wvform(:,(tet-1)*4+1 : (tet-1)*4+4, c) = mean(group(1:size(members,1),:,:));
    end

   %%% set up color code for each cluster
    linecolor = 'bgrcmyk';
    linecolor = [0 0 1; 0 1 0 ; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; .25 0 0.5; 0 0.5 0 ; 0.5 .25 0; 0.5 0 1; 0 0.5 0.5];


   %%% plot five sample waveforms for each cluster
    snipfig(tet)=figure;
   t = 1:size(X,2);
    np = 5;
    for chan = 1:4
        subplot(2,2,(chan));
        for cl = 1:n_clust;
            % for cl = 0:0
            members = find(idx==cl);
            group(1:size(members,1),1:size(X,2),1:4)=X(members,:,:);
            plot(t, group(ceil(rand(np,1)*clust_size(cl)),t,chan)'*10^6,'Color', linecolor(cl,:),'LineWidth',0.1);
            hold on
        end
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        axis([0 25 -100 50])
    end
    title(tet_title);

%%% plot histogram of spike amplitudes for all clusters
    clear hist_data;
    hist_int = 10^-6 * [-150:1:0];
    if plot_all
        histfig(tet)=figure;
        amp=zeros(n_clust,4);
        spread = amp;
        for chan = 1:4;
            subplot(2,2,chan);
            set(gcf,'DefaultAxesColorOrder',linecolor);
            for cl = 1:n_clust

                members = find(idx==cl);
                group(1:size(members,1),1:size(X,2),1:4)=X(members,:,:);
                hist_data(cl,:) = hist(min(squeeze(group(1:clust_size(cl),min_range,chan)),[],2),hist_int);
                amp(cl,chan)=mean( min(squeeze(group(1:clust_size(cl),min_range,chan)),[],2)) ;
                spread(cl,chan) = std( min(squeeze(group(1:clust_size(cl),min_range,chan)),[],2)) ;
            end

            plot(hist_int*10^6,hist_data(1:n_clust,:)','linewidth',1.5);
            axis([min(hist_int)*10^6 max(hist_int)*10^6 0 500 ]);
            %hist(min(X'),hist_int); hold on
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end
    end
    title(tet_title);
    noise
    amp
    spread

    %%% scatterplot of ICA components for all clusters (ICA1vs ICA2, ICA1 vs ICA3, etc
    clustfig(tet)=figure;
    colorstyle = ['b.g.r.c.m.y.k.'];
    n_pairs = (n_pca-1)*(n_pca)/2;
    p = 0;
    np =5000;
    N = size(score,1);
    used=zeros(size(idx));
    used(ceil(rand(np,1)*N))=1;
    for pair1 = 1:n_pca-1;
        for pair2 = pair1+1:n_pca;
            p = p+1;
            if n_pca<=4
                subplot(2,3,p);
            elseif n_pca<=6
                subplot(3,5,p)
            elseif n_pca <=8
                subplot(4,7,p);
            elseif n_pca<=12
                subplot(6,11,p)
            end

            for cl = 1:n_clust
                plot(10^6*score(find(idx==cl & used),pair1),10^6*score(find(idx==cl & used),pair2),'.','Color',linecolor(cl,:),'MarkerSize',1);
                hold on
            end
            % axis([-50 50 -50 50])
            axis off
        end
    end
    %title(tet_title);

    
    %%% calculate Lratio for each cluster
    for c= 1:n_clust
        c_score = score(find(idx==c),:);
        mn = mean(c_score,1);
        score_cov = cov(c_score);
        nonc_score = score(find(idx~=c),:);
        cov_inv = inv(cov(c_score));

        Dnon = zeros(1,size(nonc_score,1));
        Dc = zeros(1,size(c_score,1));
        var_inv = 1./var(c_score);

        for i = 1:size(c_score,1)
            Dc(i) = (var_inv.*(c_score(i,:)-mn))*(c_score(i,:)-mn)';
        end
        for i = 1:size(nonc_score,1)
            Dnon(i) = (var_inv.*(nonc_score(i,:)-mn))*(nonc_score(i,:)-mn)';
        end

        L(tet,c) = sum(1-chi2cdf(Dnon,size(score,2)));
        Lratio(tet,c) = sum(1-chi2cdf(Dnon,size(score,2)))/size(c_score,1);
        csize(tet,c) = size(c_score,1);
    end


    for c = 1:n_clust %%% waveform parameters
        c
        %      figure
        %      plot(1:size(mean_wvform,1),mean_wvform (:,(tet-1)*4+1 : (tet-1)*4+4,c)');
        title(sprintf('cluster %d',c));
        min_t = min(mean_wvform (1:15,(tet-1)*4+1 : (tet-1)*4+4,c),[],1);

        [trough_allc trig_chan] = min(min_t);
        trig_chan = (tet-1)*4+trig_chan;
        [trough(tet,c) trough_t] = min(mean_wvform (:,trig_chan,c));
        size(mean_wvform,1)
        trig_chan
        c
        trough_t
        if trough_t>=size(mean_wvform,1)
            trough_t=size(mean_wvform,1)-1;
        end
        mean_wvform ((trough_t+1):size(mean_wvform,1),trig_chan,c)
        max(mean_wvform ((trough_t+1):size(mean_wvform,1),trig_chan,c))
        [peak(tet,c) peak_t] = max(mean_wvform ((trough_t+1):size(mean_wvform,1),trig_chan,c));

        if trough_t+peak_t<size(mean_wvform,1)
            y1 = mean_wvform(trough_t+peak_t-1,trig_chan,c);
            y2 = mean_wvform(trough_t+peak_t,trig_chan,c);
            y3 = mean_wvform(trough_t+peak_t+1,trig_chan,c);
            %%% fit parabola through three points to find trough
            peak_t = peak_t - 0.5*((y3-y1)/(y3+y1-2*y2));
        end

        if trough_t>1
            y1 = mean_wvform(trough_t-1,trig_chan,c);
            y2 = mean_wvform(trough_t,trig_chan,c);
            y3 = mean_wvform(trough_t+1,trig_chan,c);
            %%% fit parabola through three points to find trough
            trough_t =  - 0.5*((y3-y1)/(y3+y1-2*y2));
        end

        halfheight = trough(tet,c)/2;
        trough_points = find(mean_wvform(:,trig_chan,c)<halfheight);
        half_minus = min(trough_points);
        if half_minus>1
            half_minus = half_minus - (halfheight - mean_wvform(half_minus,trig_chan,c))/(mean_wvform(half_minus-1,trig_chan,c) - mean_wvform(half_minus,trig_chan,c));
        end
        half_plus = max(trough_points);
        if half_plus<size(mean_wvform,1)
            halfplus = half_plus + (halfheight - mean_wvform(half_plus,trig_chan,c))/(mean_wvform(half_plus+1,trig_chan,c) - mean_wvform(half_plus,trig_chan,c));
        end
        half_width(tet,c) = half_plus - half_minus;
        peak_trough(tet,c) = peak_t-trough_t;

    end
    clear Dc Dnon nonc_score s score used
end %tet

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');
Tank_Name

clear X Xica Xold Xraw used nonc_score s score c_score

%%% save everything
[bname output_path] = uiputfile('','data folder');
bname
output_path

if bname~=0
    fname = fullfile(output_path,sprintf('cluster_data_%s_%s',Tank_Name,bname))
    save(fname);
    for t=use_tets
        fname = fullfile(output_path,sprintf('hist%s_%s_t%d',Tank_Name,bname,t));
        saveas(histfig(t),fname,'fig');
        fname = fullfile(output_path,sprintf('snip%s_%s_t%d',Tank_Name,bname,t));
        saveas(snipfig(t),fname,'fig');
        fname = fullfile(output_path,sprintf('clust%s_%s_t%d',Tank_Name,bname,t));
        saveas(clustfig(t),fname,'fig');
    end
end

toc
