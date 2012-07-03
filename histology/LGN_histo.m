function [anatomy sections] = LGN_histo(file,PaxPos);

fontSize=12;
% 
% %file names and Paxinos positions
% file = {'C:\data\histology data\ANATOMY2\dmn015_P2_A54.tif' ...
%     'C:\data\histology data\ANATOMY2\dmn015_P3_A53.tif' ...
%     'C:\data\histology data\ANATOMY2\dmn015_P4_A50.tif' ...
%     'C:\data\histology data\ANATOMY2\dmn016_P5_A53.tif'...
%     'C:\data\histology data\ANATOMY2\dmn016_P6_A53.tif'...
%     'C:\data\histology data\ANATOMY2\dmn017_P1_A55.tif'...
%     'C:\data\histology data\ANATOMY2\dmn017_P3_A49.tif'...
%     'C:\data\histology data\ANATOMY2\dmn017_P4_A53.tif'...
%     'C:\data\histology data\ANATOMY2\dmn017_P7_A52.tif'...
%     'C:\data\histology data\ANATOMY2\dmn018_P3_A50.tif'...
%     'C:\data\histology data\ANATOMY2\dmn018_P4_A53.tif'...
%     'C:\data\histology data\ANATOMY2\dmn018_P5_A54.tif'...
%     'C:\data\histology data\ANATOMY2\DMN020P1F50.tif'...
%     'C:\data\histology data\ANATOMY2\dmn021P7F50.tif'...
%     'C:\data\histology data\ANATOMY2\dmn022P10F48.tif'...
%     'C:\data\histology data\ANATOMY2\dmn023P5F53.tif'...
%     'C:\data\histology data\ANATOMY2\DMN024P9F53.tif'...
%     'C:\data\histology data\ANATOMY2\DMN024P10F51.tif'...
%     'C:\data\histology data\ANATOMY2\DMN025P5F53.tif'...
%     'C:\data\histology data\ANATOMY2\DMN025P7F52.tif'};
%     
%   PaxPos=[54 53 50 53 53 55 49 53 52 50 53 54 50 50 48 53 53 51 53 52];
%   AllenPos=[82 81 77 81 81 83 76 81 80 77 81 82 77 77 75 81 81 79 81 80];
%call drawAtlas function to get paxinos figures and 3d rendering 
[atlasFig Render3D Xall Yall Zall sections] = drawAtlas;

%read in all files 
for i = 1:length(file);
    rgbImage=imread(file{i});
    
%Find AP position. Position AP is approximated using Paxinos Atlas.
Paxinos=PaxPos(i);
%Allen=AllenPos(i);

PaxList = [44 45 46 47 48 49 50 51 52 53 54 55];
PaxAP = [-1.58 -1.70 -1.82 -1.94 -2.06 -2.18 -2.30 -2.46 -2.54 -2.70 -2.8 -2.92];
AllenList=[73 74 75 76 77 78 79 80 81 82 83];
AllenAP=[-1.855 -1.955 -2.055 -2.155 -2.255 -2.355 -2.48 -2.555 -2.78 -2.88 -2.98];
    
index = find(PaxList==Paxinos);
AP = PaxAP(index);   
paxSection=index;

% Display the original image.
composite=figure;
subplot(3, 4, 1);
imshow(rgbImage);
drawnow; % Make it display immediately showing approximate position relative to Bregma.
title(AP, 'FontSize', fontSize);


% Extract out the color bands from the original image
redBand = rgbImage(:, :, 1);
blueBand = rgbImage(:, :, 3);

lgnTrace=(blueBand==255 & redBand<255);
electrodeTrace=(redBand==255 & blueBand<255);

%plot the lgn trace
subplot(3,4,2);
imshow(~lgnTrace);
title('LGN trace', 'FontSize', fontSize);

%plot the electrode trace
subplot(3,4,3);
imshow(~electrodeTrace);
title('Electrode trace', 'FontSize', fontSize);


% x,y coordinates of lgn trace into a matrix and calculate centroid
uMperPix=1;
[YCoord XCoord] = find(flipud(lgnTrace));
CentroidX = mean(XCoord);
CentroidY = mean(YCoord);

%plot centroid
subplot(3,4,4);
plot( XCoord,YCoord);
hold on;
plot (CentroidX, CentroidY, 'cx');
title('Centroid', 'FontSize', fontSize);
axis equal;


%find tips of electrode and coordinates of electrode trace into a matrix
[YElectrode, XElectrode] = find(flipud (electrodeTrace));

[minY index]=min(YElectrode);
minX=XElectrode(index);

[maxY index]=max(YElectrode);
maxX=XElectrode(index);


% find length and angle of electrode
theta = atan2((maxY-minY),(maxX-minX));
electrodeLength=(maxY-minY)/sin(theta);
electrodeSpacing=25;

%find sites along electrode
for site=1:32;
    x(site)=minX+cos(theta)*electrodeSpacing*site;
    y(site)=minY+sin(theta)*electrodeSpacing*site;
end
 
x=fliplr(x);
y=fliplr(y);

% % rotate image so that electrode is vertical
[XCoord YCoord]=historotate(XCoord, YCoord,CentroidX, CentroidY, theta);
[XElectrode YElectrode]=historotate(XElectrode, YElectrode,CentroidX, CentroidY, theta);
[x y]=historotate(x, y,CentroidX, CentroidY, theta);
CentroidX=0;
CentroidY=0;

%plot traces with centroid and sites.
figure;
plot(XCoord, YCoord,'b.');
hold on;
plot(x, y, 'ro');
plot (XElectrode, YElectrode, 'g.');
plot (CentroidX, CentroidY, 'cx');
title(sprintf('file %d AP %f',i,AP));
axis([-600 600 -600 600]);
axis square;



%calculate site positions relative to centroid and define matrices.
Centroid = [CentroidX; CentroidY];
ElectrodeSite = [x; y];
LGN=[YCoord, XCoord];
Electrode=[YElectrode, XElectrode];


%anatomy data
anatomy(i).AP = AP;
anatomy(i).siteXY = ElectrodeSite;
anatomy(i).LGN=LGN;
anatomy(i).Electrode=Electrode;
anatomy(i).Centroid=Centroid;
anatomy(i).section = paxSection;

%plot LGN traces vs. Paxinos traces
index=find(Zall==AP*1000);
figure(composite);
subplot(3, 4, 5);
plot(XCoord, YCoord,'b.');
hold on;
plot(Xall(index), Yall(index), 'g.');


% %plot on 3d figure with electrodes in different colors
figure(Render3D);
hold on;
color = {'r.','g.','c.','m.','y.','k.'};
plot3(ones(size(x))*AP*1000,x, y, color{mod(i,6)+1});
plot3(AP*1000,x(end), y(end),'k*')
%axis([-600 600 -600 600 -3000 -1500])
axis equal
% xlabel('medial/lateral');
% ylabel('dorsal/ventral');
% zlabel('anterior/posterior');
axis off

end

%rotate 3D
for az = -30 :10:330;
    view(az,0);
    axis equal
    drawnow;
    pause(0.1)
    
end

% get analysis filename
% save(analyisfile,'anatomy','-append');


