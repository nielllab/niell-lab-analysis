    
%[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
%if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end

pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan = input('number of chans : ');
movement = input('movement data 0/1 : ');

tic
if movement==0
    flags =struct('lfpTseries',1,'lfpSpectra',1);
else
    flags =struct('lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
end


tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);

        [fname, pname] = uigetfile('*.mat','cluster data');
        load(fullfile(pname,fname));
        for i =1:length(Block_Name);
            sprintf('%d : %s ',i,Block_Name{i})
        end
        block = input('which block to analyze ? ');
        Block_Name = Block_Name{block}
        [afname, apname] = uigetfile('*.mat','analysis data');
        noisepname = apname;
        afile = fullfile(apname,afname);
        load(afile);
        use_afile=1;
        cells;

cell_range = 1:size(cells,1)

for cell_n =1:cell_range;
    ch= cell_n;
    channel_no = cells(cell_n,1)
    
    %%%
    lfp = tdtData.spectData{channel_no};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
    
    H = fspecial('average',[4 4])
    lfpnorm = imfilter(lfpnorm,H);
    
    figure
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]);
    axis xy
    df = median(diff(tdtData.spectF{channel_no}));
    dt = median(diff(tdtData.spectT{channel_no}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    
       
    if movement
        hold on
        tsamp = tdtData.mouseT;
        vsmooth = tdtData.mouseV;
        
        %plot(tsamp,(vsmooth/1.3-40),'g');
        plot(tsamp,(vsmooth/.2-40),'g');
        axis([0 max(tsamp) -40 80/df])
        
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpsc',psfilename,'-append');
    end
    title(sprintf('channel = %d',channel_no));
    
   
    %%%%

    theta = mean(lfpnorm(:,ceil(7/df):ceil(10/df)),2);
    gamma = mean(lfpnorm(:,ceil(55/df):ceil(65/df)),2);
    
    v_interp = interp1(tsamp,vsmooth,tdtData.spectT{channel_no});

    %     figure
    %     plot(v_interp,gamma(t),'o');
    %     figure
    %     plot(v_interp,theta(t),'o');
    Smean = mean(lfpnorm,2)';
    stationary = find(v_interp<2 & Smean<(5*median(Smean)));
    moving = find(v_interp>2.05  & Smean<(5*median(Smean)));

    figure
    plot(mean(lfpnorm(stationary,:),1));
    hold on
    plot(mean(lfpnorm(moving,:),1),'g');
    axis([0 70/df 0 1.2*max(mean(lfpnorm(moving,:)))]);
    set(gca,'XTick',(10:10:70)/df);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','70'})
    
    title(sprintf('site %d',channel_no));
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');

    close all
    %%%%
    
    freqs{channel_no} = tdtData.spectF{channel_no};
    df = median(diff(tdtData.spectF{channel_no}));
    specdata = tdtData.spectData{channel_no};
    normalizer = 1:size(specdata,2);
    normalizer = repmat(normalizer,size(specdata,1),1);
    specdata = specdata.*normalizer;
    
    tdata = tdtData.spectT{channel_no};
    
   
    for mv = 1:2
        
        mouse_resamp = interp1(tdtData.mouseT,tdtData.mouseV,tdata);
        if mv ==1
            timepts = mouse_resamp<1;
        else
            timepts = mouse_resamp>1;
        end
        mv_lfp(channel_no,mv,:)= nanmean(specdata(timepts,:));
    end
    
  
    plot(tdtData.spectF{channel_no}, squeeze(mv_lfp(channel_no,:,:)));
    ylim([0 1.5*prctile(mv_lfp(channel_no,1,:),95)])
    xlim([0 90])
    
    
    
    LFP_movement(channel_no).freqs = freqs{channel_no};
    LFP_movement(channel_no).mv_lfp = squeeze(mv_lfp(channel_no,:,:));
   
    close all
    
end


 save(afile,'LFP_movement','-append');
 
 
 
 