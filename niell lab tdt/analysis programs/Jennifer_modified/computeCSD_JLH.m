%function LFPtoCSD
% Code to compute CSD from LFP in response to 2-phase stimuli
% i.e. reversing checkerboard, grating, or full-field
% cmn 06-06

% connect to the tank and read the block
clear all;
pack
clear all;  % something in matlab memory management requires 2 clears ...???
pname = uigetdir('D:\','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nchan = input('number of channels : ');
duration = 0.5;
flags = struct('lfpTseries',1,'visStim',1,'lfpSpectra',1);
data = getTDTdata(Tank_Name,Block_Name,1:nchan,flags);
dt = median(diff(data.lfpT{1}))
phi = zeros(nchan,floor(duration/dt));

for ch = 1:nchan
    ch
    lfp = data.lfpData{ch};
    t = data.lfpT{ch};
    for ep = 1:length(data.stimEpocs)-1;
        onIndex = min(find(t>data.stimEpocs(2,ep)));
        phi(ch,:) = phi(ch,:)+lfp(onIndex:onIndex+floor(duration/dt-1));
    end;
end

csdfig_1=figure; set(gcf,'Name','shank 1');
csdfig_2=figure; set(gcf,'Name','shank 2');

%%% temporal filtering of LFP (done in frequency domain)
Sample_Interval=dt;
maxFreq = 200;
plot_duration = 0.2;  %%%
duration_samps = round(plot_duration/dt)
freq = fft(phi');
nyq = 1/(2*Sample_Interval)
freqInt = nyq/(0.5*size(phi,2))
maxfreqIndex = round(maxFreq/freqInt);
freq(maxfreqIndex:size(phi,2)-maxfreqIndex,:)=0;

figure(csdfig_1)
subplot(2,2,3)
plot(abs(freq(1:200,1:32)));

figure(csdfig_2)
subplot(2,2,3)
plot(abs(freq(1:200,33:64)));

phi = real(ifft(freq)');

% figure
% plot((phi(1:nchan,1:duration_samps)'));


LFPfig=figure;
minV = min(min(phi(:,1:duration_samps)));
maxV= max(max(phi(:,1:duration_samps)));
for c= 0:(nchan/4-1)
    if nchan==16
        subplot(2,2,c+1);
    elseif nchan==32
        subplot(4,2,c+1);
        else
            subplot(4,4,c+1);
    end
    plot((phi((4*c+1):(4*c+4),1:duration_samps)'));
    axis([0 duration_samps minV*1.1 maxV*1.1])
    
%     %plot LFP of all channels
%     figure
%     plot((phi(1:32,1:duration_samps)'));
%     axis([0 duration_samps minV*1.1 maxV*1.1])
    
    if c==0
        legend({'1','2','3','4'})
    end
end
title(Tank_Name);



meanLFP = ones(nchan,1)*mean(phi);
figure(csdfig_1)
subplot(2,2,1)
%plot((phi-meanLFP)');
plot((1:duration_samps)*dt*1000,phi(1:32,1:duration_samps)');
xlabel('msecs')
title(Tank_Name);

figure(csdfig_2)
subplot(2,2,1)
%plot((phi-meanLFP)');
plot((1:duration_samps)*dt*1000,phi(33:64,1:duration_samps)');
xlabel('msecs')
title(Tank_Name);


%%% calculate CSD with 2-site spacing (50M spaced probe, or 1-site for
%%% 25uMspace probe)

nchan=32
CSD_1 = (phi(1:nchan-4,:) + phi(5:nchan,:) - 2*phi(3:nchan-2,:))/(-4);
newc = imresize(CSD_1,[600 dt*1000*size(CSD_1,2)],'bilinear');  %%% interpolate at 10x sampling
figure(csdfig_1)
subplot(2,2,4);
imagesc(newc(:,(1:duration_samps*dt*1000)),1.5*[-prctile(CSD_1(:),95) prctile(CSD_1(:),95)]);
title(Tank_Name);


%second shank
nchan=64
CSD_2 = (phi(33:nchan-4,:) + phi(37:nchan,:) - 2*phi(35:nchan-2,:))/(-4);
newc = imresize(CSD_2,[600 dt*1000*size(CSD_2,2)],'bilinear');  %%% interpolate at 10x sampling
figure(csdfig_2)
subplot(2,2,4);
imagesc(newc(:,(1:duration_samps*dt*1000)),1.5*[-prctile(CSD_2(:),95) prctile(CSD_2(:),95)]);
title(Tank_Name);


%%%% calculate CSD with 4-site spacing (25uMspacing )
%first shank
CSD_1 = (phi(1:nchan-8,:) + phi(9:nchan,:) - 2*phi(5:nchan-4,:))/(-16);
newc = imresize(CSD_1,[600 dt*1000*size(CSD_1,2)],'bilinear');  %%% interpolate at 10x sampling
figure(csdfig_1)
subplot(2,2,2);
imagesc(newc(:,(1:duration_samps*dt*1000)),1.5*[-prctile(CSD_1(:),95) prctile(CSD_1(:),95)]);
title(Tank_Name);

%second shank
CSD_2 = (phi(33:nchan-8,:) + phi(41:nchan,:) - 2*phi(37:nchan-4,:))/(-16);
newc = imresize(CSD_2,[600 dt*1000*size(CSD_2,2)],'bilinear');  %%% interpolate at 10x sampling
figure(csdfig_2)
subplot(2,2,2);
imagesc(newc(:,(1:duration_samps*dt*1000)),1.5*[-prctile(CSD_2(:),95) prctile(CSD_2(:),95)]);
title(Tank_Name);

% 
% figure
% imagesc(newc(:,1:duration_samps*5),1.5*[-prctile(CSD_1(:),95) prctile(CSD_1(:),95)]);
% title(Tank_Name);
% 
% 
% figure
% imagesc(newc(:,1:duration_samps*5),1.5*[-prctile(CSD_2(:),95) prctile(CSD_2(:),95)]);
% title(Tank_Name);



[fname pname] =uiputfile('*.ps'); psfname=fullfile(pname,fname);
if fname~=0
    if exist(psfname,'file')==2;delete(psfname);end
    
    figure(csdfig)
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',fullfile(pname,fname),'-append');
    
    figure(LFPfig)
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',fullfile(pname,fname),'-append');
    
    ps2pdf('psfile', psfname, 'pdffile', [psfname(1:(end-2)) 'pdf']);
end

%
% resp=input('save results y/n','s')
% if resp=='y'
%     output_path=uigetdir('','data folder');
%          fname = fullfile(output_path,sprintf('CSD%s_%s',Tank_Name,Block_Name));
%         saveas(CSDfig,fname,'fig');
%         fname = fullfile(output_path,sprintf('tetLFP%s_%s',Tank_Name,Block_Name));
%         saveas(LFPfig,fname,'fig');
%
% end
