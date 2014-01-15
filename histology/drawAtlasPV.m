function [Paxinos3D Render3D Xall Yall Zall sections] = drawAtlasPV;
%function reads in paxinos traces, combines with AP position and generates
%a 3d composite of 2d slices and a 3d rendering

%finding scale 
% paxScale=imread('paxinos scale bar.tif');
% imshow(paxScale);

totalnum=5;

UmPix=1;
color = {'r.','g.','b.','c.','m.','y.','k.'}


%read in paxinos files and put coordinates in matrix

Paxinos3D=figure;
LGN = zeros(totalnum,1000,1000);

Xall = [];
Yall = [];
Zall = [];
for paxnum=1:totalnum;
    pathname='C:\data\matlab data\Wayne Matlab data\histology\PV traces';
     filename = sprintf('%d.tif',paxnum);%sprintf returns as string. lgn paxinos figures 44-55.
img=imread(fullfile(pathname,filename));

%separate blue pixels
redBand = img(:, :, 1);
blueBand = img(:, :, 3);

lgnTrace=(blueBand==255 & redBand<135);
[YCoord, XCoord] = find(flipud(lgnTrace));
YCoord=round(YCoord*UmPix);
XCoord=round(XCoord*UmPix);

%find centroid
CentroidX = round(mean(XCoord));
CentroidY = round(mean(YCoord));

%find Z coordinates

spacing=200;
AP = (paxnum-1)*spacing;

% for i=1:length(YCoord)
%     LGN(paxnum,XCoord(i)-CentroidX+500,YCoord(i)-CentroidY+500)=1;
% end

%figure 3d composite of 2d slices
figure(Paxinos3D);
hold on;
plot3(ones(size(XCoord))*AP, XCoord-CentroidX, YCoord-CentroidY,color{mod(paxnum,7)+1});
Xall = [Xall ; XCoord-CentroidX];
Yall = [Yall; YCoord-CentroidY];
Zall = [Zall; ones(size(XCoord))*AP];
sections(paxnum).coords=[ XCoord-CentroidX YCoord-CentroidY];
end

%3D rendering figure
[k v] = convhull(Zall,Xall,Yall,'simplify', true);
Render3D=figure;

% k = delaunay(Xall,Yall,Zall);

trisurf(k,Zall,Xall,Yall,'FaceColor','b','FaceAlpha', .5,'EdgeColor', 'none');
camlight; lighting phong
% figure
% isosurface(LGN)


