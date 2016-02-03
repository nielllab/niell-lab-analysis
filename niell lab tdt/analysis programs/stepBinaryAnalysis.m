close all; clear all

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
    blocknum = find(strcmp(Block_Name,blocknm));
    Block_Name = blocknm;
end
    
flags = struct('mouseOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampsb = tdtData.mouseT;
vsmoothsb = tdtData.mouseV;

figure
plot(tsampsb,vsmoothsb)

psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end %%% 

frame = frameEpocs{blocknum};
frameT = frame(2,:);
frameNum = frame(1,:);
frameRate = median(diff(frameT));

clear R
close all
binFrames=20;
for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
   % histbins = 5:10:max(tsampDark);
    frameR = histc(sp,frameT)/frameRate;
    
    for i = 1:(600/binFrames);
        cycAvg(i,c) = mean(frameR(mod(floor(frameNum/binFrames),600/binFrames)+1 ==i));
    end
    figure
    plot(cycAvg(:,c));
    
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end


figure
subplot(2,1,1)
imagesc(cycAvg',[0 22]); colormap jet; colorbar; axis xy; ylabel 'cell #'; xlabel 'time (s)'

subplot(2,1,2)
plot(cycAvg); ylim ([0 65]);ylabel 'cell #'; xlabel 'time (s)'
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

%append to analysisfile
save(afile,'cycAvg','tsampsb','vsmoothsb','frameR','-append');

%save pdf   
[f p] = uiputfile('*.pdf','pdf name');
ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);

%save block data
[f p] = uiputfile('*.mat','save block data?')
save(fullfile(p,f),'frameR','cycAvg','tsampsb','vsmoothsb');
