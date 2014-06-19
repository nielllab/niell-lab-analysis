    
%[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
%if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file
clear all

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


spectT_1= field2array_false(tdtData,'spectT');
spectF_1=field2array_false(tdtData,'spectF');
spectData_1=field2array_false(tdtData,'spectData');
lfpData_1=field2array_false(tdtData,'lfpData');

spectT_ch = spectT_1{1}(cell_range);
spectF_ch=spectF_1{1}(cell_range);
spectData_ch=spectData_1{1}(cell_range);
lfpData_ch=lfpData_1{1}(cell_range);

  %keyboard


for channel_no = 1:length(cell_range);
    
    
    %%%
    lfp = spectData_ch{channel_no};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
    
    H = fspecial('average',[3 3]);
    lfpnorm = imfilter(lfpnorm,H);
    
    lfpnorm_filt=lfpnorm(:,1:500);
    %lfpnorm_filt(:,320:330)=0.66*lfpnorm_filt(:,320:330);

%     moveLFP_bar(:,:,60:61)=0.6*moveLFP_bar(:,:,60:61);
%     moveLFP_bar=moveLFP_bar(:,:,1:80);
    
    figure
    imagesc(lfpnorm_filt',[0 prctile(lfpnorm(:),98)]);
    axis xy
    df = median(diff(spectF_ch{channel_no}));
    dt = median(diff(spectT_ch{channel_no}));
    set(gca,'YTick',(10:10:70)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70'})
    
       
    if movement
        hold on
        tsamp = tdtData.mouseT;
        vsmooth = tdtData.mouseV;
        
        %plot(tsamp,(vsmooth/1.3-40),'g');
        plot(tsamp,(vsmooth/.2-40),'g');
        axis([0 max(tsamp) -40 70/df])
        
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpsc',psfilename,'-append');
    end
    title(sprintf('channel = %d',channel_no));
    
   
    %%%%

    theta = mean(lfpnorm(:,ceil(7/df):ceil(10/df)),2);
    gamma = mean(lfpnorm(:,ceil(55/df):ceil(59/df)),2);
    
    v_interp = interp1(tsamp,vsmooth,spectT_ch{channel_no});

    %     figure
    %     plot(v_interp,gamma(t),'o');
    %     figure
    %     plot(v_interp,theta(t),'o');
    Smean = mean(lfp,2)';
    stationary = find(v_interp<2 & Smean<(5*median(Smean)));
    moving = find(v_interp>2.05  & Smean<(5*median(Smean)));

    
    
    figure
    plot(mean(lfpnorm_filt(moving,:),1),'g');hold on
    plot(mean(lfpnorm_filt(stationary,:),1));
    hold on
    axis([0 60/df 0 1.8*max(mean(lfpnorm_filt(moving,:)))]);
    set(gca,'XTick',(10:10:70)/df);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70'})
    
    
    title(sprintf('site %d',channel_no));
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');

    close all
    %%%%
    
 freqs{channel_no} = spectF_ch{channel_no};
   % df = median(diff(spectF_ch{channel_no}));
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
        mv_lfp_med(channel_no,mv,:)= nanmedian(specdata(timepts,:));
        
    end
   figure
   plot(tdtData.spectF{channel_no}, squeeze(mv_lfp(channel_no,:,:)));
   axis([0 70 0 1.1*max(mean(lfpnorm(moving,:)))]);
   
   figure
   plot(tdtData.spectF{channel_no}, squeeze(mv_lfp(channel_no,:,:)));
   axis([0 70 0 1.1*max(mean(lfpnorm(moving,:)))]);
    
   figure
   plot(tdtData.spectF{channel_no}, squeeze(mv_lfp_med(channel_no,:,:)));
   axis([0 70 0 1.1*max(mean(lfpnorm(moving,:)))]);
    
    LFP_movement(channel_no).freqs = freqs{channel_no};
    LFP_movement(channel_no).mv_lfp = squeeze(mv_lfp(channel_no,:,:));
    LFP_movement(channel_no).mv_lfp_med = squeeze(mv_lfp_med(channel_no,:,:));
    LFP_movement(channel_no).specdata = specdata(channel_no,:,:);
    

   
    
    
    figure
    plot(LFP_movement(channel_no).freqs,LFP_movement(channel_no).mv_lfp)
    x=squeeze(mv_lfp(1,:,1:80));
    
    figure
    imagesc(LFP_movement(channel_no).specdata);
    axis xy
   
    close all
    
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(afile,'LFP_movement','-append');
 

 
 