
function eyes = eye_detection_move(cfile,block,Block_Name, Tank_Name,redo)
close all
dbstop if error

load(cfile)%,'eyes','Block_Name','Tank_Name');

% block = find(strcmp(block,Block_Name));
% Block_Name = Block_Name(block)


if ~exist('eyes','var') | length(eyes)<block | isempty(eyes(block)) | redo
   % try
        load(cfile);  %%% contains eye images in variable 'data'
        
        data = squeeze(data);
        warning off;
        
        flags = struct('mouseOn',1,'cameraOn',1,'visStim',1);
        tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
        tsamp = tdtData.mouseT;
        vsmooth = tdtData.mouseV;
        frames = tdtData.frameEpocs;
       
    
        % to align tdt TTL with cam frames
        cameraTTL = tdtData.cameraT;
        if length(cameraTTL)>size(data,3)   %%% trim off extra TTLs after camera has been stopped
            cameraTTL = cameraTTL(1:size(data,3));
        end
        
        figure
        plot(diff(cameraTTL));
        xlabel('frame #'); ylabel('secs');
        ylim([ 0 0.2])
        
        thresh = 0.85; %pupil threshold for binarization
        puprange = [8 55]; %set range of pupil radius
        pupercent = 0.80; %set range pupil radius window
        pupchange = 0.3; %acceptable percent change in radius per framerange
        framerange = 10; %number of frames to smooth over
       % user input to select center and right points
        sprintf('Please select pupil center and top, eyeball top and right points, darkest part of eyeball')
        h1 = figure('units','normalized','outerposition',[0 0 1 1])
        imshow(data(:,:,1000))
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
        clear rad
        %convert from uint8 into doubles and threshold, then binarize       
        tic
        centroid = nan(size(data,3),2);
        rad = nan(size(data,3),1);
        centroid(1,:) = [horiz vert];
        rad(1,1) = puprad;
%%%% to test timing of stim trials%%%%
       % dbstop
        rad = squeeze(mean(mean(ddata,2),1));
        figure
        plot(rad)
%         
        
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
        
        if n>framerange && (isnan(rad(n-1)) | isnan(rad(n))) %if it's a nan or preceeded by all nans don't change puprange
            puprange = puprange;
        elseif n>framerange && (abs(1 - rad(n)/nanmean(rad(n-framerange:n-1))) > pupchange) %if % change is bigger than specified don't change puprange
            puprange = puprange;
        elseif n>framerange && (rad(n)>nanmean(rad(n-framerange:n-1))) %if radius goes up, shift range up
            puprange = puprange + round(rad(n) - nanmean(rad(n-1)));
        elseif n>framerange && (rad(n)<nanmean(rad(n-framerange:n-1))) %if radius goes down, shift range down
            puprange = puprange - round(nanmean(rad(n-1)) - rad(n));
        else
            puprange = puprange;
        end
        toc
        
        %%% plot results
        figure
        plot(cameraTTL,rad,'-b');
        hold on
        plot(cameraTTL,centroid(:,1),'g.');
        plot(cameraTTL,centroid(:,2),'r.');
        legend('radius','x','y');
        
        startT = max(tsamp(1),cameraTTL(1))
        endT = min(tsamp(end),cameraTTL(end))
        
        dt = 0.5;
        pts = startT:dt:endT;
        
        vInterp = interp1(tsamp,vsmooth,pts);
        try
            rInterp = interp1(cameraTTL(1:length(rad')),rad,pts); %if cameraTTL is longer than rad (acquisition stopped first)
        catch
            rInterp = interp1(cameraTTL,rad(1:length(cameraTTL)),pts); %if rad is longer
        end

    
     frameTime = frames(2,:);
     frameNum = frames(1,:);
   
     %interpolating frames and radius to get rad at each frame for duration of movie
     % % try
     fInterpR = interp1(cameraTTL,rad,frameTime); %if cameraTTL is longer than rad (if acquisition stopped before tdt stopped)
     fInterpX = interp1(cameraTTL,centroid(:,1),frameTime); 
     fInterpY = interp1(cameraTTL,centroid(:,2),frameTime);
     fInterpV = interp1(tsamp,vsmooth,frameTime);

%      dbstop
figure
plot(frameTime,fInterpR); xlabel('secs');
set(gcf,'Name', 'frame time & rad')

figure
plot(frameNum,fInterpR,'.'); xlabel('frame #');
set(gcf,'Name', 'frame # & rad')

for f = 1:75
    rAvg(f) = nanmean(fInterpR(mod(frameNum,75)==f-1));
end

figure
plot(rAvg);title('cyc avg pre rise');

figure
plot(frameTime(1:1000),fInterpR(1:1000))

clear R
for tr = 1:floor(max(frameNum)/75)
    for t = 1:180
        R(tr,t) = nanmean(fInterpR(frameNum == (tr-1)*75 + t));%
    end
%     R(tr,:) = R(tr,:) - R(tr,1);
end
%dbstop

figure
imagesc(R)

clear tr t
for tr = 1:floor(max(frameNum)/75)
    for t = 1:180
        V(tr,t) = nanmean(fInterpV(frameNum == (tr-1)*75 + t));%
    end
%     V(tr,:) = V(tr,:) - V(tr,1);
end
figure
imagesc(V)


%dbstop
trialSamp = fInterpR(frameNum(1:75:length(frameNum)));
figure
plot(trialSamp)

% test2= test(1:75:length(frameNum))
%   SF(1:length(lowSFpre))=1; SF(length(lowSFpre)+1:length(preSF)) =2;


%[Y,E] = discretize(frameNum,)% divides the range of X into N uniform bins, and also returns the bin edges E.


figure
plot(vInterp,'g'); hold on; plot(rInterp)
legend('velocity','radius')


        eyes(block).t = cameraTTL;
        eyes(block).x = centroid(:,1);
        eyes(block).y = centroid(:,2);
        eyes(block).rad = rad;
        eyes(block).vInterp = vInterp;
        eyes(block).rInterp = rInterp;
        eyes(block).frameT = frameTime;
        eyes(block).frameNum = frameNum;
        eyes(block).fInterpX = fInterpX;
        eyes(block).fInterpY = fInterpY;
        eyes(block).fInterpR = fInterpR;
        eyes(block).fInterpV = fInterpV;
        eyes(block).trials = R;
        eyes(block).trialV = V;
        eyes(block).trialSamp = trialSamp;
%     catch
%         
%         eyes(block).t = NaN;
%         eyes(block).x = NaN;
%         eyes(block).y = NaN;
%         eyes(block).rad = NaN;
%         eyes(block).vInterp = NaN;
%         eyes(block).rInterp = NaN;
%         
%     end
end

%%


% thresh_velocity = 1.0; %%or use 1.3, 1.5 or 2 to determine whether it changes speed distribution
% figure
% plot(tsamp,vsmooth);
% 
% %%??
% plot_duration=.5; %in second
% 
% hist_int = 0.05;
% hist_range=[0:hist_int:plot_duration];
% axis_range=[0 plot_duration 0 25];
% max_events=50000;
% 
% framerate=30
% stim_duration =.5;
% wait_duration = 2;  %duration after each stimulus
% %blank_interval = 0.1; %% length of time after stimulus not to use in calculating spontaneous rate
% %fft_int = .05;  %%% size of bins to use for fourier analyis (in sec)
% tempfreq = 0;    %%% temporal frequency of stimulus
% %blank_stim = 1; %%% is there an extra stimulus to measure spontaneous
% %full_field = 1;  %%% is there a full-field flicker?

% dbstop

% h4 = figure
% if block==1
% vidObj = VideoWriter('pre_detection.avi');
% else 
% vidObj = VideoWriter('post_detection.avi');
% 
% end
% %vidObj.FrameRate = 60;
% open(vidObj);
% 
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
%     
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
% end
% close(vidObj);
% dbstop


eyes = eyes(block)



%