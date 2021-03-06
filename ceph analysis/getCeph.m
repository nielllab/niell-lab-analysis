close all
clear all
pack
%%% load data from either .mat or .avi
[f p] = uigetfile({'*.mat;*.avi'});
if strcmp(f(end-3:end),'.mat')  %%% .mat file
    load(fullfile(p,f));
    im = squeeze(graypairtwo);
    for fr = 1:size(im,3)
        if fr==1
            imwrite(im(:,:,fr),fullfile(p,[f(1:end-3) '.tif']),'tiff');
        else
            imwrite(im(:,:,fr),fullfile(p,[f(1:end-3) '.tif']),'tiff','Writemode','append');
        end
    end    
else %%% .avi file
   display('reading avi')
   vr= VideoReader(fullfile(p,f));
    fr=0;
    im = zeros(vr.Height,vr.Width,vr.Framerate*vr.duration,'uint8');
   % im = zeros(vr.Height,vr.Width,500,'uint8');
    while (hasFrame(vr)) 
        fr=fr+1;
        frame = readFrame(vr,'native');
        im(:,:,fr) = frame(:,:,1);
    end
    
end
im = single(im);
im(im==0)=1;

figure
imagesc(im(:,:,1))
title('im')

%%% check stability
dx = 4;
imsparse = double(im(dx:dx:end,dx:dx:end,:));
imsparse = reshape(imsparse,size(imsparse,1)*size(imsparse,2),size(imsparse,3));
figure
imagesc(corrcoef(imsparse));

%%% get brightness as mean across image
brightness = mean(imsparse,1);
figure
plot(brightness)
clear imsparse

%%% normalize image for mean brightness
display('normalizing')

for fr = 1:size(im,3);
    im(:,:,fr) = im(:,:,fr)/brightness(fr);
end

figure
imagesc(im(:,:,1))
title('norm_im')

%%% mask image regions that don't change over the whole session(background)
display('masking')
std_img = std(im(:,:,dx:dx:end),[],3);
figure
imagesc(std_img);
figure
imagesc(std_img>0.05);

for fr = 1:size(im,3);
   im(:,:,fr) = im(:,:,fr).*(std_img>0.05);
end
im(im==0)=NaN;
figure
imagesc(im(:,:,1));
title('mask im')

%%% choose dark regions (this threshold may need to be adjusted
thresh_im = im<0.8;
figure
imagesc(thresh_im(:,:,1))


%%% binary operation to remove holes
display('closing holes')
se = strel('disk',5);
for fr = 1:size(thresh_im,3)
    thresh_im(:,:,fr) = imfill(imclose(thresh_im(:,:,fr),se),'holes');
end
figure
imagesc(thresh_im(:,:,1))


%%% go through binarized images and grab the biggest one
display('getting regions')
%labelmov = zeros(size(thresh_im,1),size(thresh_im,2),3,size(thresh_im,3));
clear obj
for fr = 1:size(thresh_im,3)
    label = bwlabel(thresh_im(:,:,fr),8);
    props = regionprops(label,'Area','Orientation','Centroid');
    sz = [props.Area];
    [m ind] = max(sz);
    obj(fr) = props(ind);
    orient(fr) = props(ind).Orientation;
    centroid(fr,:) = props(ind).Centroid;
    label = label==ind;
    im(:,:,fr)=label.*im(:,:,fr);
   % labelmov(:,:,1,fr) = label(:,:,fr); labelmov(:,:,2,fr) = label(:,:,fr); labelmov(:,:,3,fr) = label(:,:,fr);
end
figure
imagesc(label)


figure
plot(centroid(:,1),centroid(:,2));axis ([1 size(im,2) 1 size(im,1)])
title('centroid position')

figure
imagesc(im(:,:,1));
clear label thresh_im


%%% check bounds of centroid
w = 250;
centroid = max(centroid,w+1);
centroid(:,1) = min(centroid(:,1),size(im,2)-w);
centroid(:,2) = min(centroid(:,2),size(im,1)-w);
centroid = round(centroid);

%%% choose block around centroid and rotate it to match
display('getting blocks')
range = -w:w;
imblock_rot = zeros(length(range),length(range),size(im,3));
for fr = 1:size(im,3);
    imblock = im(centroid(fr,2)+range,centroid(fr,1)+range,fr);
    imblock_rot(:,:,fr) = imrotate(imblock,90-orient(fr),'crop');
    if fr>1;
      %%% orientation is ambiguous (up/down) so use correlation to first
      %%% image as match
      if sum(sum(flipud(imblock_rot(:,:,fr)).*imblock_rot(:,:,1)))>sum(sum((imblock_rot(:,:,fr)).*imblock_rot(:,:,1)))
            imblock_rot(:,:,fr) = flipud(imblock_rot(:,:,fr));
        end
    end
    
end

figure
imagesc(imblock_rot(:,:,1))
axis equal

%%% choose region aligned with body
imblock_rot = imblock_rot(:,100:399,:);

%%% save out aligned and rotated images
imblock_rot_sm = imresize(imblock_rot,0.25);
for fr = 1:size(imblock_rot,3)
    fr
    if fr==1
        imwrite(imblock_rot_sm(:,:,fr),fullfile(p,[f(1:end-3) '_align.tif']),'tiff');
    else
        imwrite(imblock_rot_sm(:,:,fr),fullfile(p,[f(1:end-3) '_align.tif']),'tiff','Writemode','append');
    end
end

clear im

%%% calculate correlation of image over time, to find stable periods
%%% (just for trouble-shooting, doesn't get used
dx = 2;
imsparse = imblock_rot(dx:dx:end,dx:dx:end,:);
imsparse = reshape(imsparse,size(imsparse,1)*size(imsparse,2),size(imsparse,3));
imsparse = imsparse>0;
tic
cc = corrcoef(imsparse);
toc
figure
imagesc(cc)

%%% get correlation of each frame with mean image, to find good matches
meanimg = (mean(imblock_rot>0,3)>0.5); axis equal
for fr = 1:size(imblock_rot,3)
    match(fr) = sum(sum((imblock_rot(1:250,:,fr)>0) .* meanimg(1:250,:)));  %%% 1:250 selects top half of image
end
match = match/ (sum(sum(meanimg(1:250,:))));

figure
plot(match)
figure
hist(match,0.01:0.02:1)

%%% find stationary blocks, corr>corrthresh for longer than mindur frames
mindur = 20; corrthresh  =0.97;
movement = [0 find(match<corrthresh) length(match+1)];
duration = diff(movement);
stable = find(duration>mindur);
duration(stable)
movement(stable)
movement(stable+1)
matchfilt= zeros(size(match));
for i = 1:length(stable)
    matchfilt(movement(stable(i))+2:movement(stable(i))+duration(stable(i))-2)=1;
end
figure
plot(matchfilt,'g');
hold on
plot(match)

%%% save out aligned tif movie for stable periods
imblock_rot_sm = imresize(imblock_rot(:,:,matchfilt>0.95),1);
% for fr = 1:size(imblock_rot_sm,3)
%     fr
%     if fr==1
%         imwrite(imblock_rot_sm(:,:,fr),fullfile(p,[f(1:end-4) '_stable.tif']),'tiff');
%     else
%         imwrite(imblock_rot_sm(:,:,fr),fullfile(p,[f(1:end-4) '_stable.tif']),'tiff','Writemode','append');
%     end
% end
% 
%%% "warp" moving along x-axis, aligning to mean. Intended to compensate
%%% for side to side movement/twisting
% mn = mean(imblock_rot_sm>0,3);
% display('warping')
% for fr = 1:size(imblock_rot_sm,3)
%     for i = 1:size(imblock_rot_sm,1);
%         [xc lag] = xcorr(double(mn(i,:)),double(squeeze(imblock_rot_sm(i,:,fr)>0)));
%         [m ind] = max(xc);
%         shift = lag(ind);
%         im_warp(i,:,fr) = circshift(squeeze(imblock_rot_sm(i,:,fr))',[shift 0]);
%     end
% end

%%% save out "warped" movie
im_warp = imblock_rot_sm;
for fr = 1:size(imblock_rot_sm,3)
    fr
    if fr==1
        imwrite(im_warp(:,:,fr),fullfile(p,[f(1:end-4) '_warp.tif']),'tiff');
    else
        imwrite(im_warp(:,:,fr),fullfile(p,[f(1:end-4) '_warp.tif']),'tiff','Writemode','append');
    end
end

%%% hardcoded section to analyze!!!

display('might want to check warped movie to determine frames')
keyboard

im_warp = im_warp(:,:,1:end);

%%% select chromatophores
mn = mean(im_warp,3);
mask = zeros(size(im_warp));
clear sz
maskfig=figure
mnfig= figure;
imagesc(mean(im_warp,3)); axis equal; colormap gray;
for spot = 1:10
    spot
    mnfig= figure;
    imagesc(mean(im_warp,3)); axis equal; colormap gray;
    [y x]  = ginput(2);
    chrom(spot,1)=y(1); chrom(spot,2)=x(1);
    thresh = (mn(round(x(1)),round(y(1))) + mn(round(x(2)),round(y(2))))/2;
    thresh
    
    for fr=1:size(im_warp,3)
        obj =  bwselect(im_warp(:,:,fr)<thresh,y(1),x(1),4);
       outline = bwperim(obj);
        mask(:,:,fr) =   mask(:,:,fr)+spot*outline;      
        sz(fr,spot) = sqrt(sum(sum(obj)));
    end
    close(mnfig)
    
    figure
    imagesc(mean(mask,3)); axis equal
end

%%% plot sizes
figure
plot(sz);
sz (sz==0)= NaN; sz(sz>15)=NaN;
legend; ylabel('size'); xlabel('frame')
figure
plot(sz);

%%% plot correlations
cc= corrcoef(sz,'rows','pairwise');
figure
imagesc(cc,[0 1]); colormap jet
title('size correlations')

%%% draw connectivty figure
mainfig = figure;
subplot(2,4,[1 5])
imagesc(mn); colormap(gray);axis equal
hold on;
for i = 1:size(sz,2);
    plot(chrom(i,1),chrom(i,2),'o','Linewidth',4);
    for j= 1:size(sz,2);
        if cc(i,j)>0.4;
            plot([chrom(i,1) chrom(j,1)],[chrom(i,2) chrom(j,2)],'b','Linewidth',12*(cc(i,j)-0.4),'Color',cmapVar(cc(i,j),0.4,0.8,jet));
        end
    end
end
axis([60 240 50 500]); axis off
set(gca,'LooseInset',get(gca,'TightInset'))

%%% normalize sizes and subtract min size

clear normsz
for sp = 1:size(sz,2);
    normsz(:,sp) = (sz(:,sp)-min(sz(:,sp)))/nanstd(sz(:,sp));
    normsz(:,sp) = naninterp(normsz(:,sp));
end

%%% plot normalized sizes
subplot(2,4,2);
plot(normsz)

%%% PCA analysis
[coeff score latent] = pca(normsz,'rows','pairwise');
figure(mainfig)
subplot(2,4,3)
plot(latent/sum(latent)); hold on; plot(latent/sum(latent),'o')
set(gca,'LooseInset',get(gca,'TightInset'))
figure
imagesc(coeff,[-0.75 0.75]); colormap jet

%%% non-negative matrix factorization
for n = 1:10;
    [w h err(n)] = nnmf(normsz(5:end,:),n,'replicates',50);
end
figure
plot(err)

tic
[w h] = nnmf(normsz(5:end,:),3,'replicates',50);
toc
% figure
% imagesc(w)
% figure
% imagesc(h')

%%% plot nnmf components
figure(mainfig)
loc = [6 7 8]
for c = 1:size(h,1)
 subplot(2,4,loc(c));
    hold off ; imagesc(mn); colormap(gray);axis equal; hold on
    for i = 1:size(sz,2);
        if h(c,i)>0
            plot(chrom(i,1),chrom(i,2),'g.','Markersize',50*h(c,i), 'Linewidth',12*h(c,i));
        end
    end
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end

%%% plot nnmf timecourse
subplot(2,4,4);
plot(w)


%%% save out tif movie with chromatophore outlines
for fr = 1:size(mask,3)
    fr
    im = mat2im(im_warp(:,:,fr),gray);
    im(:,:,2) = im(:,:,2)+mask(:,:,fr);
    overlay=im;
    if fr==1
        imwrite(im,fullfile(p,[f(1:end-4) '_mask.tif']),'tiff');
    else
        imwrite(im,fullfile(p,[f(1:end-4) '_mask.tif']),'tiff','Writemode','append');
    end
end



