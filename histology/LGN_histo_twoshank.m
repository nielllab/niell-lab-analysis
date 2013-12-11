function [anatomy sections] = LGN_histoPV(file,PaxPos);

fontSize=12;

%call drawAtlas function to get paxinos figures and 3d rendering
[atlasFig Render3D Xall Yall Zall sections] = drawAtlasPV;

if size(file,2)==2;
    twoshank=1;
else
    twoshank=0;
end

file
%read in all files
for i = 1:length(file);
   sprintf('filenum %d',i)
   if isempty(file{i,1})
        anatomy(i).AP = [];
        anatomy(i).siteXY = [];
        anatomy(i).LGN=[];
        anatomy(i).Electrode=[];
        anatomy(i).Centroid=[];
        anatomy(i).section = [];
    else
        for shank=1:2
            sprintf('%d shank',shank)
            if shank==1 || (shank==2 && ~isempty(file{i,2}))
                
                rgbImage=imread(file{i,shank});
                
                
                if i==6
                    keyboard
                end
                for c = 1:3
                    rgbImage(:,:,c) = fliplr(rgbImage(:,:,c));
                end
                
                
                %Find AP position. Position AP is approximated using Paxinos Atlas.
                Paxinos=PaxPos(i);
                %Allen=AllenPos(i);
                
                spacing = .200;
                PaxList = [1 2 3 4 5 ];
                PaxAP = 0:spacing:(length(PaxList)-1)*spacing;
                
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
                
                anatomy(i).LGN=LGN;
                anatomy(i).Electrode=Electrode;
                anatomy(i).Centroid=Centroid;
                anatomy(i).section = paxSection;
                if shank==1
                    anatomy(i).siteXY = ElectrodeSite;
                else
                    anatomy(i).siteXY = [anatomy(i).siteXY ElectrodeSite];
                end
                
                %plot LGN traces vs. Paxinos traces
                AP*1000
                unique(Zall)
                index=find(Zall==round(AP*1000));
                size(index)
                figure(composite);
                subplot(3, 4, 5);
                plot(XCoord, YCoord,'b.');
                hold on;
                plot(Xall(index), Yall(index), 'g.');
                axis([-500 500 -500 500])
                axis equal
                
                
                % %plot on 3d figure with electrodes in different colors
                
                figure(Render3D);
                hold on;
                color = {'r.','g.','c.','m.','y.','k.'};
                plot3(ones(size(x))*AP*1000,x, y, color{mod(i,6)+1});
                %plot3(AP*1000,x(end), y(end),'k*')
                %axis([-600 600 -600 600 -3000 -1500])
                if i==1
                    axis equal
                end
                % xlabel('medial/lateral');
                % ylabel('dorsal/ventral');
                % zlabel('anterior/posterior');
                axis off
                
                % LGNmovie(i)=getframe(gcf);
                
            end
        end
    end
end

% vid = VideoWriter('LGNmovie slow');
% vid.FrameRate=6;
% open(vid);
% writeVideo(vid,LGNmovie);
% close(vid)

%rotate 3D
for az = -30 :10:330;
    view(az,0);
    axis equal
    drawnow;
    pause(0.1)
    
end

% get analysis filename
% save(analyisfile,'anatomy','-append');


