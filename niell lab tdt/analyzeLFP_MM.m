close all
clear all

pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

% Tank_Name='06222011_mlr_stim_awake'
% Block_Name='wn3_72hz'
nChan=16;
tic
flags =struct('lfpTseries',1,'lfpSpectra',1);

tdtData= getTDTdata(Tank_Name, Block_Name, 1:16, flags);
toc

%ranges  = [5 10; 20 30; 50 60]; %%% modify this for different freqs
chans = [1 5 12 16];
for c=1:length(chans);
   ch = chans(c);
   lfp = tdtData.spectData{ch};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
   
    figure   
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]);
    axis xy
    df = median(diff(tdtData.spectF{ch}));
    dt = median(diff(tdtData.spectT{ch}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    title(sprintf('channel = %d',ch));
    
    spect = median(lfpnorm,1);
    freqs = tdtData.spectF{ch};
    figure
    plot(tdtData.spectF{ch},spect);
    for f= 1:120;
        sp(c,f) = median(spect(freqs>f-0.5 & freqs<f+0.5));
    end
    hold on
    plot(sp(c,:));
    
end

%%% sp contains mean LFP spectrum (normalized) on 4 chans as a function of
%%% true frequency (Hz)