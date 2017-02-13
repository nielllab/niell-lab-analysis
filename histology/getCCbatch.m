close all
clear all

drct = dir('*.tif')

for filenum = 1:length(drct)

fname = drct(filenum).name;

cc = double(imread(fname,1));
edu = double(imread(fname,2));

im(:,:,1) = edu;
im(:,:,2) = cc;
im(:,:,3) = cc;
figure
subplot(2,2,1)
imshow(im/255)
title(fname)

subplot(2,2,2)
edupos = edu(edu>0);
hist(edupos(:),0:255); xlim([0 255])

eduMask = edu>150;


eduMask = imopen(eduMask,strel('disk',3));
subplot(2,2,3)
imagesc(eduMask)


eduLabel = bwlabeln(eduMask);
% figure
% imagesc(eduLabel)

eduObjs = regionprops(eduLabel);
area = [eduObjs(1:end).Area];
% figure
% hist(area)

thresh = 20;
for i = 1:length(area)
    if area(i)<thresh
        eduMask(eduLabel==i)=0;
    end
end

% figure
% imagesc(eduMask)
eduLabel = bwlabeln(eduMask);
obj = regionprops(eduMask);
length(obj)

clear ccObj
for i = 1:length(obj)
    objMask = (eduLabel==i);
    outer = imdilate(objMask,strel('disk',2));
    inner = imerode(objMask,strel('disk',2));
    circ = outer>inner;
%       if i/10 ==round(i/10)
%           figure
%         imagesc(circ);
%       end
    ccObj(i) = mean(mean(cc(circ)));
end

range = [10:20:255];
h=hist(ccObj,range); h= h/sum(h);
ccfilt = imfilter(cc,fspecial('disk',6));
ccpos = ccfilt(ccfilt>0);
h2 = hist(ccpos(:),range); h2 = h2/sum(h2);

subplot(2,2,4)
bar(range,h); hold on; plot(range,h2,'g','LineWidth',2);
xlim([0 255]) ;legend('edu','random');

d = h-h2;
sig = round(sum(d(d>0))*length(ccObj));

% st = std(ccpos(:)); mn = mean(ccpos(:));
% thr = mn + 2*st;
thr = prctile(ccpos(:),95);
sig2std =  sum(ccObj>thr);


plot([thr thr],[0 0.2],'r','LineWidth',2)


% ccThresh=input('threshold : ');
% pos = false(size(cc));
% for i =1:length(obj)
%     if ccObj(i)>ccThresh
%         pos(eduLabel==i)=1;
%     end
% end
% 
% figure
% im(:,:,3)=0;
% imshow(im)
% 
% im(:,:,3) = uint8(pos*255);
% 
% figure
% imshow(im);
% 
% figure
% imshow(squeeze(im(:,:,2)))
% 
% sprintf('%d out of %d positive based on threshold',sum(ccObj>ccThresh),length(ccObj))
% 
% sprintf('%d out of %d positive based on random',sig,length(ccObj))

sprintf('%d out of %d positive based on 95 prctile',sig2std,length(ccObj))

allfiles{filenum} = fname;
nObj(filenum) = length(ccObj);
nPos(filenum) = sig2std;
threshold(filenum) = thr;

end
[f p] = uiputfile('*.xlsx');
fname = fullfile(p,f);
delete(fname)
xlswrite(fname,allfiles');
xlswrite(fname,nObj',1,'B1');
xlswrite(fname,nPos',1,'C1');





