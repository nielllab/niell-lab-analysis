%eyetracking & neural data & running data

close all
clear all
dbstop if error

[fname, pname] = uigetfile('*.mat','cluster data');
clustfile=fullfile(pname,fname);
load(clustfile);
blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');

[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
load(afile);
[pname fname] = fileparts(afile);
Block_Name = Block_Name{blocknum}
    
[cfname, cpname] = uigetfile('*.mat','eyetracking data');
cfile = fullfile(cpname,cfname);
load(cfile);
[pname fname] = fileparts(cfile)

psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end %%% 

data = squeeze(data);   
warning off;

flags = struct('mouseOn',1,'cameraOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampBar = tdtData.mouseT;
vsmoothBar = tdtData.mouseV;

% to align tdt TTL with cam frames
cameraTTL = tdtData.cameraT;
figure
plot(diff(cameraTTL));
xlabel('frame #'); ylabel('secs');
ylim([ 0 0.2])

thresh = 0.8; %pupil threshold for binarization
puprange = [8 45]; %set

%user input to select center and right points
sprintf('Please select pupil center and top, eyeball top and right points, darkest part of eyeball')
h1 = figure('units','normalized','outerposition',[0 0 1 1])
imshow(data(:,:,10))
[cent] = ginput(5);
close(h1);
yc = cent(1,2); %pupil center y val
xc = cent(1,1); %pupil center x val
horiz = (cent(4,1) - xc); %1/2 x search range
vert = (yc - cent(3,2)); %1/2 y search range
puprad = yc - cent(2,2); %initial pupil radius
% puprange = [round(puprad - puprad*pupercent) round(puprad + puprad*pupercent)]; %range of pupil sizes to search over
ddata = double(data); 
binmaxx = cent(5,1);
binmaxy = cent(5,2);

for i = 1:size(data,3)
    binmax(i) = mean(mean(mean(ddata(binmaxy-3:binmaxy+3,binmaxx-3:binmaxx+3,i))));
end
for i = 1:size(ddata,3)
    bindata(:,:,i) = (ddata(yc-vert:yc+vert,xc-horiz:xc+horiz,i)/binmax(i) > thresh);
end

%convert from uint8 into doubles and threshold, then binarize

tic
centroid = nan(size(data,3),2);
rad = nan(size(data,3),1);
centroid(1,:) = [horiz vert];
rad(1,1) = puprad;
for n = 2:size(data,3)
    [center,radii,metric] = imfindcircles(bindata(:,:,n),puprange,'Sensitivity',0.995,'ObjectPolarity','dark');
    if(isempty(center))
        centroid(n,:) = [NaN NaN]; % could not find anything...
        rad(n) = NaN;
    else
        [~,idx] = max(metric); % pick the circle with best score
        centroid(n,:) = center(idx,:);
        rad(n,:) = radii(idx);
    end
end
toc

clear R
close all
for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
    histbins = 5:5:max(tsampBar);
    R(1:length(histbins),c) = hist(sp,histbins)/10;
   figure
    plot(histbins,R(1:length(histbins),c),'-r'); 
   hold on
    plot(tsampBar,vsmoothBar,'-g');
   plot(0.1:0.1:size(data,3)/10,rad,'-b')
   plot(0.1:0.1:size(data,3)/10,(centroid(:,1)/1.1),'.m')
   plot(0.1:0.1:size(data,3)/10,(centroid(:,2)/1.1),'.c')
   ylim([0 40]);
   legend('sp/sec','cm/sec','radius','x pos','ypos')
   set(gcf, 'PaperPositionMode', 'auto');
   print('-dpsc',psfilename,'-append');
    
end

startT = max(tsampBar(1),cameraTTL(1))
endT = min(tsampBar(end),cameraTTL(end))

dt = 0.5;
pts = startT:dt:endT;
vInterp = interp1(tsampBar,vsmoothBar,pts);
rInterp = interp1(cameraTTL,rad(1:length(cameraTTL)),pts);

figure
plot(vInterp); hold on; plot(rInterp,'r')

% %plot cross correlation of velocity and radius
[corr_vrad lags] = xcorr(rInterp-mean(rInterp),vInterp-mean(vInterp),60/dt,'coeff');
figure
plot(lags*dt,corr_vrad);
hold
plot([0 0],[0 1],'g-')

set(gcf, 'PaperPositionMode', 'auto');
 print('-dpsc',psfilename,'-append')


% figure
% scatter(Vinterp, Rinterp, '.b')
% xlabel('radius');ylabel('cm/sec')
% lsline;
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpsc',psfilename,'-append')

% %for FR and rad
%scatter(vsmoothBar,R(1:length(histbins),c))
% 
%cross correlation of FR and rad (and lags)for each unit
for c = 1:length(spikeT)
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
    histbins = [pts pts(end)+dt];
    rate=hist(sp,histbins)/dt;
    rate = rate(2:end);
    R(1:length(histbins),c) = hist(sp,histbins)/dt;
    [ corr_Frad lags] = xcorr(rate-mean(rate),rInterp-mean(rInterp),60/dt,'coeff');
    [ corr_Fvel lags] = xcorr(rate-mean(rate),vInterp-mean(vInterp),60/dt,'coeff');
    
    figure
    subplot(2,2,1)
    plot(rate/max(rate)); hold on; plot(vInterp/max(vInterp),'g'); plot(rInterp/max(rInterp),'r');
    subplot(2,2,2)
    plot(lags,corr_Frad)
    title('FR and rad')
    subplot(2,2,3)
    plot(lags,corr_Fvel)
    title('FR and velocity')
    preISI = sp(2:end-1)-sp(1:end-2);
    postISI = sp(3:end) - sp(2:end-1);
    subplot(2,2,4)
    plot(preISI,postISI,'.')
    xlabel('preISI'); ylabel('postISI')
  
   
   xcorr(rate-mean(rate),rate-mean(rate))
    %
        set(gcf, 'PaperPositionMode', 'auto');
         print('-dpsc',psfilename,'-append');
end


   
%velocity and FR
% %    scatter(vsmoothBar,R(1:length(histbins),c))        
% 
% h2 = figure
% hold on
% plot(0.1:0.1:size(data,3)/10,rad,'b-')
% plot(0.1:0.1:size(data,3)/10,(centroid(:,1)/1.1),'.m')
% plot(0.1:0.1:size(data,3)/10,(centroid(:,2)/1.1),'.c')
% hold off
% legend('radius','x pos','ypos')

%plot radius vs speed


% h3=figure
% plot(rad, 

% video with tracking
% h4 = figure
% for i = 1:size(data,3)
%  
%     subplot(1,2,1)
%     imshow(data(yc-vert:yc+vert,xc-horiz:xc+horiz,i));
%     colormap gray
%     hold on
%     circle(centroid(i,1),centroid(i,2),rad(i))
%     drawnow
%     hold off
%     
%     subplot(1,2,2)
%     imshow(bindata(:,:,i));
%     colormap gray
%     hold on
%     circle(centroid(i,1),centroid(i,2),rad(i))
%     drawnow
%     hold off 
% %     mov(i) = getframe(gcf)
% %     vid = VideoWriter('predoi_tracking.avi')
% %     open(vid); writeVideo(vid,mov); close(vid)
% end


save(afile,'tsampBar','vsmoothBar','R','centroid','rad','-append');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');
    
[f p] = uiputfile('*.pdf','pdf name');
ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);
