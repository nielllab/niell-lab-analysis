% Track_Green_mice.m Tristan Ursell/Jennifer Hoy/ Emily Galleta 4-21-17
% load up "movie" file one frame at a time, images saved as series of
    % .tif files
%pathname = 'E:\Dropbox\mucsimol_V1\Trial videos for Tristans tracking\Iphone5\test_green_ears\30Hz_1080_SS75_ISO50_Tiff\'
%pathname='X:\Tracks\BHDN'
pathname = '/Users/jennifer/Desktop/green ear videos/controls/'
%cd 'F:\GreenEarTracking\iPhone\5_25_17\'

%pathname='F:\GreenEarTracking\iPhone\5_25_17\ALGSo101_PS\'
files=dir(fullfile(pathname,'*.tif'));

%load series of .tif files as "movie"
% bin=2;c=1;
% for i=1:length(files)/bin ; %crude downsample    
% 
% [img(:,:,:,i)] = imread(files(c+bin).name); %#ok<SAGROW>
%     imshow(img(:,:,:,i));
%     drawnow   
%   c=c+bin;
% end




%load full length video
for i=1:length(files) ;

[img(:,:,:,i)] = imread(files(i).name); 
%    imshow(img(:,:,:,i));
 %    drawnow   
    
end

% 
% 
% for i=1:size(img,4)%
%    imshow(img(:,:,:,i));
%   % imshow(trim(:,:,:,500));
%    drawnow   
% end

scale=22; % pixels/cm
fps=30; % could change between videos

%% insert function for color thresholding
%%
numframes = size(img, 4);

%segment video based on color space
% clear R_vec G-vec B-vec
% R_vec=img(:,:,1);
% R_vec=R_vec(:);
% %histR=hist(R_vec,1:1:255);
%  
% G_vec=img(:,:,2);
% G_vec=G_vec(:);
% %histG=hist(G_vec,1:1:255);
%  
% B_vec=img(:,:,3);
% B_vec=B_vec(:);
%histB=hist(B_vec,1:1:255);

%%
% color based segmentation using thresholded masks,

%%apply mask to whole video
%parpool

for k=1:numframes
mask1(:,:,k)=img(:,:,2,k)>(img(:,:,1,k)+10); %G>R+10
mask2(:,:,k)=img(:,:,2,k)>(img(:,:,3,k)+10); %G>B+10
mask3(:,:,k)=mask1(:,:,k).*mask2(:,:,k);%ear mask
end


for k=1:numframes
mask4(:,:,k)=(img(:,:,1,k)>(img(:,:,3,k)+10)).*(img(:,:,1,k)<155).*(~mask3(:,:,k)); %R>B+10 & R<150 & ~mask3: is the red parts of mouse
mask5(:,:,k)=(img(:,:,1,k)<120).*(img(:,:,2,k)<120).*(img(:,:,3,k)<120).*(~mask3(:,:,k)).*(~mask4(:,:,k));%body mask
%mask6(:,:,k)=(img(:,:,2,k)>(img(:,:,3,k)+3)).*(~mask2(:,:,k));  %green greater than blue by 3 and 
%mask7(:,:,k)=mask6(:,:,k).*mask4(:,:,k);%cricket mask
end

figure
imagesc(mask4(:,:,1))
% t = bwlabel(mask4(:, :, k)); %mask for the ears
% cum = regionprops(t, 'area', 'centroid','Eccentricity');
% area_vector = [cum.Area];
% shap_vec=[cum.Eccentricity];

% %show masks across whole video
% for i=50:100;%numframes
% imshow(mask4(:,:,i));hold on
% drawnow
% end
% 
% figure
% imagesc(mask4(:,:,200));
% 
%% path of body and ear by finding object centroids
clear k centroidD centroidE AreaC CentroidC CentroidE1 CentroidE2 CentroidB
CentroidB = zeros(length(numframes), 2);
CentroidE1=zeros(length(numframes), 2);
CentroidE2=zeros(length(numframes), 2);
CentroidC=zeros(length(numframes), 2);
AreaC=zeros(length(numframes), 1);


for k = 1:numframes
LE = bwlabel(mask3(:, :, k)); %mask for the ears
sE = regionprops(LE, 'area', 'centroid');
area_vector = [sE.Area];
if size(area_vector,2)>1
    [j h]=sort(area_vector,'descend');
    CentroidE1(k,:)=sE(h(1),1).Centroid;
    CentroidE2(k,:)=sE(h(2),1).Centroid;
    
    LB = bwlabel(mask5(:, :, k)); %body mask
    sB = regionprops(LB, 'area', 'centroid');
    area_vector = [sB.Area];
    [tmp, idx] = max(area_vector);
    CentroidB(k, :) = sB(idx(1)).Centroid;
    
%     LC = bwlabel(mask7(:, :, k));
%     sC = regionprops(LC, 'area', 'centroid','Eccentricity');
%     area_vector = [sC.Area];
%     [tmp, idx] = max(area_vector);
%     CentroidC(k,:)=sC(idx(1)).Centroid;
%     AreaC(k)=sC(idx(1)).Area;
else
    CentroidE1(k,:)=NaN;
    CentroidE2(k,:)=NaN;
    CentroidB(k,:)=NaN;
    %CentroidC(k,:)=NaN;
end
end

%% get cricket track
type=1% type=1 is live type=0 equals virtual
for k = 1:numframes
   if type==1 
clear E_vec area_vector sc Lc

Lc = bwlabel(mask4(:,:,k)); %mask for the cricket
sc = regionprops(Lc, 'area', 'centroid','Eccentricity');
area_vector = [sc.Area];
%E_vec=[sc.Eccentricity];
big=find(area_vector>200 );


%determine if there are cricket sized objects in the frame
if ~isempty(big) & length(big)==1
CentroidC(k,:)=[sc(big).Centroid]
elseif ~isempty(big) & length(big)>1; %if there are more than 2 objects that meet criteria than select the cricket one manually
    sprintf 'click on centroid of cricket'
    imshow(img(:,:,:,k));
     [X,Y]=getpts;
    CentroidC(k,:)=[X(1),Y(1)];
  
else
    CentroidC(k,:)=NaN

end 
   else
       imshow(img(:,:,:,1));
     [X,Y]=getpts;
    CentroidC(1:end,1)= [X(1)];
    CentroidC(1:end,2)= [Y(1)];
 
   end
   
end

% for i=1:numframes
% imshow(img(:,:,:,i));hold on
% plot(CentroidC(i,1),CentroidC(i,2),'c*');hold on
% %plot(Right(i,1),Right(i,2),'r*');hold on
% %plot(CentroidB(i,1),CentroidB(i,2),'y*');hold on
% drawnow
% 
% %mov(i)=getframe(gcf);
% end

% figure
% imagesc(mask7(:,:,k))
% %

% plot how matlab defines the centroids for the objects of interest
figure
plot(CentroidE1(:,1),CentroidE1(:,2),'m');hold on
plot(CentroidE2(:,1),CentroidE2(:,2),'g');hold on
plot(CentroidB(:,1),CentroidB(:,2),'k');hold on
plot(CentroidC(:,1),CentroidC(:,2),'r*');hold on

%% assign each object to as either the left or right ear
clear Left Right
[Left, Right]=LeftEar(img,mask3,CentroidB,CentroidE1,CentroidE2,0); %1 in last entry position implies automatic assignment of Left ear vs. right, if "0" not auto 1st click is left ear, then Shift+ click the right ear center

%check for correct assignment of ear ID should print out automatically when
%you run LeftEar

% figure
%     plot(Left(:,1),Left(:,2),'g');hold on
%     plot(Right(:,1),Right(:,2),'r');hold on
%    %current
%     plot(Left(end,1),Left(end,2),'k*');hold on
%     plot(Right(end,1),Right(end,2),'r*');hold on
%     %previous
%     plot(Left(end-1,1),Left(end-1,2),'g*');hold on
%     plot(Right(end-1,1),Right(end-1,2),'m*');hold on

% for i=900:1000%numframes
% imshow(img(:,:,:,i));hold on
% %plot(Left(i,1),Left(i,2),'c*');hold on
% %plot(Right(i,1),Right(i,2),'r*');hold on
% plot(CentroidC(i,1),CentroidC(i,2),'r*');hold on
% drawnow
% %mov(i)=getframe(gcf);
% end


%% determine cricket track and exclude frames where mouse is almost in
% contact with cricket

%find frames where mouse touch cricket

for i=1:numframes
    DistLE_C(i)=sqrt((Right(i,1) - CentroidC(i,1)).^2 + (Right(i,2) - CentroidC(i,2)).^2); 
    DistB_C(i)=sqrt((CentroidB(i,1) - CentroidC(i,1)).^2 + (CentroidB(i,2) - CentroidC(i,2)).^2); 
end

% touch= find(AreaC<150);

touch= find(DistLE_C<50 |DistB_C<50 );
j=2
figure
imagesc(img(:,:,:,j)); hold on
plot(CentroidC(j,1),CentroidC(j,2),'c*');hold on

%omit cricket from frames where mouse is touching it
CentroidCNT(:,:)=CentroidC(:,:);
CentroidCNT(touch,:)= NaN;

figure
plot(Left(:,1),Left(:,2),'r');hold on
plot(Right(:,1),Right(:,2),'g');hold on
plot(CentroidB(:,1),CentroidB(:,2),'c'); hold on
plot(CentroidCNT(:,1),CentroidCNT(:,2),'k');

%fill in blank cricket data with previous location of cricket

CentroidCNT = fill_nans(CentroidCNT);
% Replaces the nans in each column with 
% previous non-nan values.

% figure
% plot(Left(:,1),Left(:,2),'r');hold on
% plot(Right(:,1),Right(:,2),'g');hold on
% plot(CentroidB(:,1),CentroidB(:,2),'c'); hold on
% plot(CentroidC(:,1),CentroidC(:,2),'k');
%% generate movie of tracked points

for i=1:numframes
imshow(img(:,:,:,i));hold on
plot(Left(i,1),Left(i,2),'g*');hold on
plot(Right(i,1),Right(i,2),'y*');hold on
plot(CentroidB(i,1),CentroidB(i,2),'r*');hold on
plot(CentroidCNT(i,1),CentroidCNT(i,2),'c*');hold on

drawnow
%mov(i)=getframe(gcf);
end

% vidObj = VideoWriter('tracked.avi');
%         vidObj.FrameRate = 10;
%         open(vidObj);
%         writeVideo(vidObj,mov);
%         close(vidObj);




% make video 

% for i=1:numframes
% imshow(mask3(:,:,i));hold on
% plot(Left(i,1),Left(i,2),'g*');hold on
% plot(Right(i,1),Right(i,2),'m*');hold on
% plot(CentroidB(i,1),CentroidB(i,2),'y*');hold on
% plot(CentroidCNT(i,1),CentroidCNT(i,2),'r*');hold on
% drawnow
% 
% mov(i)=getframe(gcf);
% 
% end
% 
% vidObj = VideoWriter('tracked.avi');
%         vidObj.FrameRate = 10;
%         open(vidObj);
%         writeVideo(vidObj,mov);
%         close(vidObj);

%% plot final data set as a "movie"
for i=1:numframes;
imshow(img(:,:,:,i));hold on
plot(Left(i,1),Left(i,2),'g*');hold on
plot(Right(i,1),Right(i,2),'m*');hold on
plot(CentroidCNT(i,1),CentroidCNT(i,2),'r*');hold on
drawnow
if i==numframes;
plot(Left(:,1),Left(:,2),'g');hold on
plot(Right(:,1),Right(:,2),'m');hold on
plot(CentroidCNT(:,1),CentroidCNT(:,2),'r');
end

% mov(i)=getframe(gcf);

end
% 
% vidObj = VideoWriter('tracked.avi');
%         vidObj.FrameRate = 10;
%         open(vidObj);
%         writeVideo(vidObj,mov);
%         close(vidObj);


%% build in manual input interface for color thresholding
%%
%color space of pixels, validate seperate of colorspace on images
% figure;
% subplot(2,2,1)
% hold on
% rv1=G_vec>(R_vec+10);
% plot(R_vec(rv1),G_vec(rv1),'.','color',[0.7 0.7 0])
% plot(R_vec(~rv1),G_vec(~rv1),'x','color',[1 1 0])
% plot([0 255],[0 255],'k--')
% xlabel('R')
% ylabel('G')
% axis equal
% xlim([0 255])
% ylim([0 255])
% 
% subplot(2,2,2)
% hold on
% plot(R_vec,B_vec,'.','color',[0.7 0 0.7])
% plot([0 255],[0 255],'k--')
% xlabel('R')
% ylabel('B')
% axis equal
% xlim([0 255])
% ylim([0 255])
% %box on
% 
% subplot(2,2,3)
% hold on
% bv1=G_vec>(B_vec+10);
% plot(G_vec(bv1),B_vec(bv1),'.','color',[0 0.7 0.7])
% plot(G_vec(~bv1),B_vec(~bv1),'.','color',[0 1 1])
% plot([0 255],[0 255],'k--')
% xlabel('G')
% ylabel('B')
% axis equal
% xlim([0 255])
% ylim([0 255])

% %check that the correct objects are segmented based on color space split
% subplot(2,2,4)
% test2=uint8(sqrt(double(img(:,:,1)).^2+double(img(:,:,2)).^2+double(img(:,:,3)).^2));
% test3(:,:,1)=test2.*uint8(~mask3(:,:,1)).*uint8(~mask5(:,:,1));
% test3(:,:,2)=test2.*uint8(~mask4(:,:,1)).*uint8(~mask5(:,:,1));
% test3(:,:,3)=test2.*uint8(~mask3(:,:,1)).*uint8(~mask4(:,:,1));
% imagesc(test3)

%% Goal 4:  write functions that use the 'annotated' ears to calculate ear
%centroids (regionprops), head centroid. generate and save out track data
 

%%

%cd 'F:\GreenEarTracking\iPhone\5_25_17\F:\GreenEarTracking\iPhone\5_25_17\ALGSo101_PS'
% filename = 'TracksFile.mat';
% save(filename,'Left','Right','CentroidB', 'CentroidC');
%         
% pathname='F:\GreenEarTracking\iPhone\5_25_17\ALGSo101_PS\'

[Tfname Tpname] = uiputfile('*.mat','Tracks file');
fullaname = fullfile(Tpname, Tfname)
save(fullaname, 'Left','Right','CentroidB', 'CentroidCNT','touch')


