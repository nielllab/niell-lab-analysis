clear all; close all

if ~exist('clustfile','var')  %%%stand alone run
[fname, pname] = uigetfile('*.mat','cluster data');
clustfile=fullfile(pname,fname);
load(clustfile);
blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
  load(afile);
    [pname fname] = fileparts(afile);
    Block_Name = Block_Name{blocknum}
else   %%% if using batch
    load(clustfile)
    load(afile)
    blocknm = find(strcmp(Block_Name,blocknm));
    Block_Name = blocknm;
end

flags = struct('mouseOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampDark = tdtData.mouseT;
vsmoothDark = tdtData.mouseV;

dt=0.5

figure
plot(tsampDark,vsmoothDark)

psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end %%% 


clear R
close all
% for c = 1:length(spikeT)
%     sp = spikeT{c};
%     sp = sp-(blocknum-1)*10^5;
%     sp = sp(sp>0 & sp<10^5);
%     histbins = 5:10:max(tsampDark);
%     darkR(1:length(histbins),c) = hist(sp,histbins)/10;
%   % [corr_Fvel lags] = xcorr(darkR{c}-mean(darkR{c}),vsmoothDark-mean(vsmoothDark),60/dt,'coeff')
%     figure
%     subplot(2,1,1)

%     plot(histbins,darkR(1:length(histbins),c));
%     hold on
%     plot(tsampDark,vsmoothDark,'g');
%     ylim([0 15]);
%     legend('sp/sec','cm/sec')
%     %sp = spikeT{c};
%     subplot(2,1,2)
%     hist(spikeT{c},0:10^5:max(spikeT{c}))
%     %subplot(2,2,3)
%     %plot(lags,corr_Fvel)
%     %title('FR and velocity')
%     
% %     set(gcf, 'PaperPositionMode', 'auto');
% %     print('-dpsc',psfilename,'-append');
%     blockSpike{c} = sp;
%
% end

for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
    histbins = 5:10:max(tsampDark);
    rate=hist(sp,histbins)/dt
    rate = rate(2:end);
    
    darkR(1:length(histbins),c) = hist(sp,histbins)/dt;
    % darkR(:,c) = hist(sp,histbins)/10;
    figure
    subplot(1,2,1)
    % plot(rate/max(rate)); hold on; plot(vInterp/max(vInterp),'g'); plot(rInterp/max(rInterp),'r');
    
    plot(histbins,darkR(1:length(histbins),c))    
    hold on
    plot(tsampDark,vsmoothDark,'g');
    ylim([0 20]);
    legend('sp/sec','cm/sec')
    hold off
    subplot(1,2,2)
    
    [corr_Fvel lags] = xcorr(rate-mean(rate),vsmoothDark-mean(vsmoothDark),60/dt);
    plot(lags,corr_Fvel)
    title('FR and velocity')
    % title_text = sprintf('ch%d c%d',channel_no,clust_no);
    
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    blockSpike{c} = sp;
end


save(afile,'tsampDark','vsmoothDark','darkR','-append');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');
    
[f p] = uiputfile('*.pdf','pdf name');
ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);

[f p] = uiputfile('*.mat','save block data?')
save(fullfile(p,f),'blockSpike','tsampDark','vsmoothDark', 'darkR','corr_Fvel');
