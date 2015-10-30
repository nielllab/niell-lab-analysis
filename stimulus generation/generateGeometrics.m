triBig = imread('geometric shapes-05.jpg');
triSml = imread('geometric shapes-40.jpg');
sqrBig = imread('geometric shapes-15.jpg');
sqrSml = imread('geometric shapes-30.jpg');
crcBig = imread('geometric shapes-10.jpg');
crcSml = imread('geometric shapes-35.jpg');
starBig = imread('geometric shapes-20.jpg');
starSml = imread('geometric shapes-25.jpg');

shapes = zeros(256,256,8); %8 basic starting shapes, 4 of 2 sizes each
shapes(:,:,1) = triBig;
shapes(:,:,2) = triSml;
shapes(:,:,3) = sqrBig;
shapes(:,:,4) = sqrSml;
shapes(:,:,5) = crcBig;
shapes(:,:,6) = crcSml;
shapes(:,:,7) = starBig;
shapes(:,:,8) = starSml;

stimuli = 127*ones(256,256,128,144);
asp = size(stimuli,3);
leg = zeros(144,3);

cnt = 1;
for a = 1:size(shapes,3)
    for color = 1:2
        if color==1
            stim = shapes(:,:,a);
            stim(stim<=127) = 127; stim(stim>127) = 0;
            for direc = 1:9
                if direc==1
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = stim;
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==2;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==3;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2,2);
                    end        
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==4;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2,1);
                    end 
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==5;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==6;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2-n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==7;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2+n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==8;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2+n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==9;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2-n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                end
            end
        else 
            stim = shapes(:,:,a);
            stim(stim<=127) = 127; stim(stim>127) = 255;
            for direc = 1:9
                if direc==1
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = stim;
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==2;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==3;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2,2);
                    end        
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==4;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2,1);
                    end 
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==5;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==6;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2-n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==7;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2+n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                 elseif direc==8;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,-asp/2+n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==9;
                    for n = 1:asp
                        stimuli(:,:,n,cnt) = circshift(stim,asp/2-n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                end
            end
        end
    end
end

                    
% figure
% for j = 1:128
%     imshow(stimuli(:,:,j,6));
% end