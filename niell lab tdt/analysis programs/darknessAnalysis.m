close all
clear all

[fname, pname] = uigetfile('*.mat','cluster data');
clustfile=fullfile(pname,fname);
load(clustfile);
blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
  load(afile);
    [pname fname] = fileparts(afile);
    Block_Name = Block_Name{blocknum}
    
flags = struct('mouseOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampDark = tdtData.mouseT;
vsmoothDark = tdtData.mouseV;

figure
plot(tsampDark,vsmoothDark)

psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end %%% 

clear R
close all
for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
    histbins = 5:10:max(tsampDark);
    darkR(:,c) = hist(sp,histbins)/10;
    figure
    plot(histbins,darkR(:,c));
    hold on
    plot(tsampDark,vsmoothDark,'g');
    ylim([0 15]);
    legend('sp/sec','cm/sec')
    
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    
    
end

figure
plot(histbins,mean(darkR,2));hold on
plot(tsampDark,vsmoothDark,'g');
ylim([0 15]);
legend('sp/sec','cm/sec'); title('average')

save(afile,'tsampDark','vsmoothDark','darkR','-append');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');
    
[f p] = uiputfile('*.pdf','pdf name');
ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);