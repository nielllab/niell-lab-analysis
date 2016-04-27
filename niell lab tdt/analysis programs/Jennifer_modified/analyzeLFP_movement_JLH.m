%function analyzeLFP_movement
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block

[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


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

for ch = 1:nChan;
   
    lfp = tdtData.spectData{ch};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
    
    H = fspecial('average',[6 6])
    lfpnorm = imfilter(lfpnorm,H);
    
    figure
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]); xlim ([0 400]); ylim ([0 500]);
    axis xy;
    df = median(diff(tdtData.spectF{ch}));
    dt = median(diff(tdtData.spectT{ch}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'}); 
    xlabel 'time (s)' ; ylabel 'Frequency (Hz)';
    
    if movement
        hold on
        tsamp = tdtData.mouseT;
        vsmooth = tdtData.mouseV;
        
        %plot(tsamp,(vsmooth/1.3-40),'g');
        plot(tsamp,(vsmooth/1-40),'g');
        axis([0 max(tsamp) -40 80/df])
        title(sprintf('channel = %d',ch));
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpsc',psfilename,'-append');    
    end
    
    %%%%
    
    theta = mean(lfpnorm(:,ceil(7/df):ceil(10/df)),2);
    gamma = mean(lfpnorm(:,ceil(55/df):ceil(65/df)),2);
    
    v_interp = interp1(tsamp,vsmooth,tdtData.spectT{ch});
    
    %     figure
    %     plot(v_interp,gamma(t),'o');
    %     figure
    %     plot(v_interp,theta(t),'o');
    Smean = mean(lfpnorm,2)';
    stationary = find(v_interp<0.3 & Smean<(5*median(Smean)));
    moving = find(v_interp>0.35  & Smean<(5*median(Smean)));
    
    figure
    plot(mean(lfpnorm(stationary,:),1));
    hold on
    plot(mean(lfpnorm(moving,:),1),'g');
    axis([0 70/df 0 1.2*max(mean(lfpnorm(moving,:)))]);
    set(gca,'XTick',(10:10:80)/df); 
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'});
    xlabel 'Frequency (Hz)';
    title(sprintf('channel = %d',ch));
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

close all


ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);


