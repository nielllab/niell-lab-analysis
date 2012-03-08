tet_title = 'tetrode'
global goodcells;
    plot_all = 1;
    plot_chan =0;
   
threshold_voltage=0;

for tet=1

    clear X
%%%% Bass's data
%     X= shiftdim(WV,2)/1000;   
 %   X = X(:,25:49,:);
  %idx_all = zeros(16,size(WV,3),'uint8');

%%%% Willie's data

% for i = 1:5
%      X(:,:,i) = reshape(squeeze(x(:,i)),41,size(x,1)/41)';
%  end
%  for i = 1:4
%     event_times_all(i,1:size(X,1))=X(:,1,1);
%  end
%  figure
%  plot(X(1:30,:,2)');
%  idx_all = zeros(16,size(X,1),'uint8');
%  X= X(:,14:38,2:5);


load next4
idx_all = zeros(16,size(X,1),'uint8');
X=X/(1.8*10^6);
event_times_all = event_times_all*(.05*10^-3);


    clear Xshift Xold etimes_old newtimes;
    pack
    Xsize= size(X);


    use_pca=1;
    individual_score=0;



    N=size(X,1)
    used = 1:N;
    trigT = 6;
    shiftrange=6;

    sz = 25-shiftrange;
    Xshift = zeros(N,sz,4);

    
    %%% remove snippets that have an initial voltage that is beyond the
    %%% usual noise range, as they are generally artifacts
 
    v0thresh = 80*10^-6;    %%% threshold for initial point
    vmax_thresh = 800*10^-6;  %%% threshold for largest excursion

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

    for chan = 1:4
        % [Y I] = min(squeeze(X(:,trigT:trigT+shiftrange,chan)),[],2);
        [Y I] = min(squeeze(mean(X(:,trigT:trigT+shiftrange,:),3)),[],2);
        for i = 1:N
            Xshift(i,:,chan) = X(i,I(i) : (I(i)+sz-1),chan);

        end
    end
    figure
    hist(I);
    X= Xshift(used,:,:);

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

    choose_pca=1;

    if choose_pca
        n_pca=8;
    else
        n_pca =8;
    end
    mn = mean(X,3);
    pc_data = [X(:,:,1) X(:,:,2) X(:,:,3) X(:,:,4)];
    size(pc_data)
    [s coeff u] = fastica(pc_data'*10^6,'numOfIC',n_pca,'lastEig',n_pca,'stabilization','on','g','tanh','approach','symm');
    score = s';
    clear s;
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

    clear pc_data Xraw I mn
    pack
    clear pca_use
    figure(ica_fig)

    goodcells = zeros(1,n_pca);

    if choose_pca
        %     ginput(1)
        %     n_pca=input('how many ica to use?  ');
        %     for i = 1:n_pca;
        %         pca_use(i) = input('enter next ica : ')
        %     end
        done=0;
        pause
       pca_use = find(goodcells);
        n_pca = length(pca_use);
%         while ~done
%             k = waitforbuttonpress;
%             if k==1 %%% keyboard
%                 done=1;
%             end
%         end


    else
        pca_use=1:n_pca;
    end
    %pca_use=1:n_pca
    score = score(:,pca_use);

    pca_use
    fname='pca.fet.1'
    dlmwrite(fname,n_pca);
    dlmwrite(fname,score,'-append','delimiter',' ');
    sprintf('clustering data .....')
    tic; dos('C:\KlustaKwik_amd pca 1 -MinClusters 5 -MaxClusters 10 -MaxPossibleClusters 12 -Verbose 0 -nStarts 2','-echo'); toc
    %n_pca=4;
    cl_name = 'pca.clu.1';
    idx = uint8(dlmread(cl_name));
    n_clust = idx(1)
    if sum(idx==1)<10   %%% if no members of outlier cluster
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


    for c = 1:n_clust
        members = find(idx==c);
        group(1:size(members,1),1:size(X,2),1:4)=X(members,:,:);

        mean_wvform(:,(tet-1)*4+1 : (tet-1)*4+4, c) = mean(group(1:size(members,1),:,:));
    end

    snipfig(tet)=figure;
    linecolor = 'bgrcmyk';
    linecolor = [0 0 1; 0 1 0 ; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; .25 0 0.5; 0 0.5 0 ; 0.5 .25 0; 0.5 0 1; 0 0.5 0.5];

    %t = Sample_Interval*(0:size(group1,2)-1);
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




    clear hist_data;
    hist_int = 10^-6 * [-250:1:0];
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


    amp
    spread

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




clear X Xica Xold Xraw used nonc_score s score c_score

[bname output_path] = uiputfile('','data folder');
use_tets=1;
Tank_Name = 'Willie'
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
