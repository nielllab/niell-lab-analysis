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
for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
   % histbins = 5:10:max(tsampDark);
    frameR = histc(sp,frameT)/frameRate;
    
    for i = 1:600;
        cycAvg(i) = mean(frameR(mod(frameNum,600)+1 ==i));
    end

    figure
    plot(cycAvg);


% cellspikes = frameR(i,:,:);
% normR(c,:,:) = cellspikes/max(cellspikes(:));
end

%allR = [squeeze(normR(:,:,1)) squeeze(normR(:,:,2))];
figure
imagesc(cycAvg(i));

dur=600
figure 
imagesc(cycAvg);
% axis xy;
% hold on
% plot(tsampsb/dt,5*vsmoothsb/max(vsmoothsb)-5 ,'g'); ylim([-5 size(frameR,1)+0.5]);
% %plot spike times for each unit
% for c = 1:length(spikeT)
%     sp = spikeT{c};
%     figure
%     hist(spikeT{c})
% end
% 
% figure
% plot(histbins,mean(darkR(1:length(histbins),2)));
% hold on
% plot(tsampDark,vsmoothDark,'g');
% ylim([0 15]);
% legend('sp/sec','cm/sec'); title('average')



save(afile,'tsampDark','vsmoothDark','darkR','-append');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');
    
[f p] = uiputfile('*.pdf','pdf name');
ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);

[f p] = uiputfile('*.mat','save block data?')
save(fullfile(p,f),'blockSpike','tsampDark','vsmoothDark');
