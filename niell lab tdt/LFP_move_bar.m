    
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
movement = 1;

tic

flags =struct('lfpTseries',1,'lfpSpectra',1,'mouseOn',1);

[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
load(afile);
use_afile=1;

if exist('LFP_movement','var')
 clear LFP_movement 
end


cells;

cell_range = 1:size(cells,1)
cell_range=(cells(:,1));

tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);


spectT_1= field2array(tdtData,'spectT');
spectF_1=field2array(tdtData,'spectF');
spectData_1=field2array(tdtData,'spectData');
lfpData_1=field2array(tdtData,'lfpData');

spectT_ch = spectT_1{1}(cell_range);
spectF_ch=spectF_1{1}(cell_range);
spectData_ch=spectData_1{1}(cell_range);
lfpData_ch=lfpData_1{1}(cell_range);

  


for channel_no = 1:length(cell_range);
    
    
    %%%
    lfp = spectData_ch{channel_no};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
    
    H = fspecial('average',[4 4])
    lfpnorm = imfilter(lfpnorm,H);
    
    figure
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]);
    axis xy
    df = median(diff(spectF_ch{channel_no}));
    dt = median(diff(spectT_ch{channel_no}));
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
    
    v_interp = interp1(tsamp,vsmooth,spectT_ch{channel_no});

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
    
    freqs{channel_no} = spectF_ch{channel_no};
    df = median(diff(spectF_ch{channel_no}));
    specdata = spectData_ch{channel_no};
    normalizer = 1:size(specdata,2);
    normalizer = repmat(normalizer,size(specdata,1),1);
    specdata = specdata.*normalizer;
    
    tdata = spectT_ch{channel_no};
    
   
    for mv = 1:2
        
        mouse_resamp = interp1(tdtData.mouseT,tdtData.mouseV,tdata);
        if mv ==1
            timepts = mouse_resamp<1;
        else
            timepts = mouse_resamp>1;
        end
        mv_lfp(channel_no,mv,:)= nanmean(specdata(timepts,:));
    end
    
  
    plot(spectF_ch{channel_no}, squeeze(mv_lfp(channel_no,:,:)));
    ylim([0 1.5*prctile(mv_lfp(channel_no,1,:),95)])
    xlim([0 90])
    
    
    
    LFP_movement(channel_no).freqs = freqs{channel_no};
    LFP_movement(channel_no).mv_lfp = squeeze(mv_lfp(channel_no,:,:));
   
    close all
    
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(afile,'LFP_movement','-append');
 

 
 