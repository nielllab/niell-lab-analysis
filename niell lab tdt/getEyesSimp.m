function eyes = getEyesSimp(clustfile,afile,cfile,block,redo)
%%% Extract eye tracking data for a given block
%%% Checks to see if data already exists
%%% If not, reads in camera data and asks for points
%%% Results are saved to analysis file in variable 'eyes', with an entry for each block
%%% Returns eye data for the selected block

load(clustfile,'Block_Name','Tank_Name');
%blocknum = find(strncmp(block,Block_Name,3));
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}


if ~exist('eyes','var') | length(eyes)<blocknum | isempty(eyes(blocknum)) | redo
    try
        load(cfile);  %%% contains eye images in variable 'data'
        
        data = squeeze(data);
        warning off;
        
        flags = struct('cameraOn',1);
       % tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
        
        % to align tdt TTL with cam frames
       % cameraTTL = tdtData.cameraT;
       % figure
       % plot(diff(cameraTTL));
        %xlabel('frame #'); ylabel('secs');
        %ylim([ 0 0.2])
        
        thresh = 0.8; %pupil threshold for binarization
        puprange = [2 25]; %set range of pupil radius
        
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
        
figure;imagesc(squeeze(bindata(:,:,150)));
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
        plot(rad,'b.');
        hold on
        plot(centroid(:,1),'g.');
        plot(centroid(:,2),'r.');
        legend('radius','x','y');
        
      %  eyes(blocknum).t = cameraTTL;
        eyes(blocknum).x = centroid(:,1);
        eyes(blocknum).y = centroid(:,2);
        eyes(blocknum).rad = rad;
        
    catch
        
       % eyes(blocknum).t = NaN;
        eyes(blocknum).x = NaN;
        eyes(blocknum).y = NaN;
        eyes(blocknum).rad = NaN;
    end
    
    save(cfile,'eyes','-append');
    
end

eyes = eyes(blocknum);  %%% return data for this specific block