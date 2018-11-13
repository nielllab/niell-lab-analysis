%%
clear all

cd('C:\Users\nlab\Box Sync\Phil Niell Lab\Behavior\Stimuli')

xsz = 1280; %set to 1x mag on psych stim controller
ysz = xsz*9/16;
dist = 20;
width = 50;
widthdeg = 2*atand(0.5*width/dist);
degperpix = widthdeg/xsz;

rad1 = [5]; %radius of circles in degrees
rad2 = [10]; %radius1 and 2 must be same length

xpos1 = 0.3; xpos2 = 0.7; ypos = 0.5;
rad1x = round(xpos1*xsz); rad2x = round(xpos2*xsz);
rady = round(ypos*ysz);

stim1 = 128*ones(ysz,xsz);stim2=stim1;

[X,Y] = meshgrid(1:xsz,1:ysz);

dist = sqrt((X-rad1x).^2 + (Y-rady).^2);
stim1(dist<(rad1/degperpix)) = 0;
dist = sqrt((X-rad2x).^2 + (Y-rady).^2);
stim1(dist<(rad2/degperpix)) = 0;
figure;imshow(stim1)

dist = sqrt((X-rad1x).^2 + (Y-rady).^2);
stim2(dist<(rad2/degperpix)) = 0;
dist = sqrt((X-rad2x).^2 + (Y-rady).^2);
stim2(dist<(rad1/degperpix)) = 0;
figure;imshow(stim2)

imwrite(stim1,'stim1.bmp','BMP')
imwrite(stim2,'stim2.bmp','BMP')

save 2AFCgrating
