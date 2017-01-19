function eyes = getEyes_angie(cfile, block,redo)
%%% Extract eye tracking data for a given block
%%% Checks to see if data already exists
%%% If not, reads in camera data and asks for points
%%% Results are saved to camera file in variable 'eyes', with an entry for each block
%%% Returns eye data for the selected block

load(cfile,'eyes','Block_Name','Tank_Name');

blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name(blocknum)

%load(afile,'eyes');

if ~exist('eyes','var') | length(eyes)<blocknum | isempty(eyes(blocknum)) | redo
    try
        load(cfile);  %%% contains eye images in variable 'data'
        
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
        puprange = [10 55]; %set range of pupil radius
        
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
        
        %%% plot results
        figure
        plot(rad,'-b');
        hold on
        plot(centroid(:,1),'g.');
        plot(centroid(:,2),'r.');
        legend('radius','x','y');
        
        startT = max(tsampBar(1),cameraTTL(1))
        endT = min(tsampBar(end),cameraTTL(end))
        
        dt = 0.5;
        pts = startT:dt:endT;
        
        vInterp = interp1(tsampBar,vsmoothBar,pts);
%         try
        %rInterp = interp1(cameraTTL(1:length(rad')),rad,pts); %if cameraTTL is longer than rad (due to blinks?)
%         catch
        rInterp = interp1(cameraTTL,rad(1:length(cameraTTL)),pts); %if rad is longer
%         end
        
       
        figure
        plot(vInterp,'g'); hold on; plot(rInterp)
        legend('velocity','radius')
    

        eyes(blocknum).t = cameraTTL;
        eyes(blocknum).x = centroid(:,1);
        eyes(blocknum).y = centroid(:,2);
        eyes(blocknum).rad = rad;
        eyes(blocknum).vInterp = vInterp;
        eyes(blocknum).rInterp = rInterp;
  
    catch
        
%         eyes(blocknum).t = NaN;
%         eyes(blocknum).x = NaN;
%         eyes(blocknum).y = NaN;
%         eyes(blocknum).rad = NaN;
%         eyes(blocknum).vInterp = NaN;
%         eyes(blocknum).rInterp = NaN;
%         
    end
end
 save(cfile,'eyes','-append');
eyes = eyes(blocknum);  %%% return data for this specific block
