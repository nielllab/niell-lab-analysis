function [Left Right]=LeftEar(img,Cbody,Cear1,Cear2,auto) 
dbstop if error
if auto ==1

numframes=length(Cbody);
CE1_left=zeros(length(Cbody),1);
for i=1:numframes-1
    
    if Cbody(i,2)<Cbody(i+1,2)-1;
        mousemov=1; %1=down
    elseif Cbody(i,2)>Cbody(i+1,2)+1;
        mousemov=2; %2=up
    elseif Cbody(i,1)< Cbody(i+1,1)-1
        mousemov=3; %3=right, no up or down motion
    elseif Cbody(i,1)> Cbody(i+1,1)+1;
        mousemov=4;%4=left, no up or down motion
    else
        mousemov=5;%mouse is stationary
    end
    
    if mousemov==1 && Cear1(i,1)-Cear2(i,1)<0;
        CE1_left(i)=0;
    elseif mousemov==1 && Cear1(i,1)-Cear2(i,1)>0;
        CE1_left(i)=1;
    elseif mousemov==2 && Cear1(i,1)-Cear2(i,1)<0;
        CE1_left(i)=0;
    elseif mousemov==2 && Cear1(i,1)-Cear2(i,1)>0;
        CE1_left(i)=1;
    elseif mousemov==3 && Cear1(i,1)-Cear2(i,1)<0;
        CE1_left(i)=0;
    elseif mousemov==3 && Cear1(i,1)-Cear2(i,1)>0;
        CE1_left(i)=1;
        
    elseif mousemov==4 && Cear1(i,1)-Cear2(i,1)<0;
        CE1_left(i)=1;
    elseif mousemov==4 && Cear1(i,1)-Cear2(i,1)>0;
        CE1_left(i)=0;
    elseif mousemov==5 && i~=1
        CE1_left(i)=CE1_left(i-1);
    else
        CE1_left=NaN;
        
    end
    
end
        CE1_left(end+1)=CE1_left(end-1);      
        
       Left=zeros(numframes, 2);
       Right=zeros(numframes, 2);
             
       for i=1:numframes
           if i==1
               sprintf 'click left ear then hold shift and click right ear'
               imshow(img(:,:,:,i));
               [X,Y]=getpts;
               Left(1,:)=[X(1),Y(1)];
               Right(1,:)=[X(2),Y(2)];
         
           elseif i>1 & ~isnan(CE1_left(i))
               if CE1_left(i)==1
                   Left(i,:)=Cear1(i,:);
                   Right(i,:)=Cear2(i,:);
               else
                   Left(i,:)=Cear2(i,:);
                   Right(i,:)=Cear1(i,:);
               end
           elseif i>1 & isnan(CE1_left(i));
               imshow(img(:,:,:,i));
               sprintf 'click left ear then hold shift and click right ear'
               [X,Y]=getpts;
               Left(i,:)=[X(1),Y(1)];
               Right(i,:)=[X(2),Y(2)];
              
           end
       end

% for i=1:numframes;
% imshow(mask(:,:,i+1));hold on
% plot(Left(i+1,1),Left(i+1,2),'g*');hold on
% plot(Right(i+1,1),Right(i+1,2),'y*');hold on
% plot(Cbody(i+1,1),Cbody(i+1,2),'c*');hold on
% drawnow
% %mov(i)=getframe(gcf);
% end

%% check that Left(i) X,Y is closer to Left(i-1) X,Y than Right(i-1)X,Y
    for i=2:numframes;
        clear DistLL and DistLR
    %dist between  Left(i) and Left(i-1); r = sqrt((head(:,1) - targ(:,1)).^2 + (head(:,2) - targ(:,2)).^2)
    DistLL = sqrt((Left(i,1) - Left(i-1,1)).^2 + (Left(i,2) - Left(i-1,2)).^2);
    DistLR = sqrt((Left(i,1) - Right(i-1,1)).^2 + (Left(i,2) - Right(i-1,2)).^2);
   
    if DistLL > DistLR +10 & Left(i,:)==Cear1(i,:);
       Left(i,:)=Cear2(i,:);
       Right(i,:)=Cear1(i,:);
    elseif DistLL > DistLR +10 & Left(i,:)==Cear2(i,:);
       Left(i,:)=Cear1(i,:);
       Right(i,:)=Cear2(i,:);
    elseif DistLL > DistLR | DistLL+10 > DistLR;
        sprintf 'click left ear then hold shift and click right ear'
        imshow(img(:,:,:,i));
        [X,Y]=getpts;
        Left(i,:)=[X(1),Y(1)];
        Right(i,:)=[X(2),Y(2)];
    end
    
    end

%   figure
%   for j=1:numframes
%       plot(Left(j,1),Left(j,2),'r*');hold on
%       plot(Right(j,1),Right(j,2),'g*');hold on
%     drawnow
%     keyboard
%   end
    
figure
plot(Left(:,1),Left(:,2),'r');hold on
plot(Right(:,1),Right(:,2),'g');hold on
plot(Cbody(:,1),Cbody(:,2),'c')

% find 1st point for left and right ears manually, left ear is the first
% selected point, right ear is the second selected point, hold shift when
% selecting 2nd pt to store them out.
else
    numframes=length(Cbody);
    Left= Cear1;
    Right=Cear2; 
    
    sprintf 'click left ear then hold shift and click right ear'
    imshow(img(:,:,:,1));
    [X,Y]=getpts;
    
    Left(1,:)=[X(1),Y(1)];
    Right(1,:)=[X(2),Y(2)];

%check that Left(i) X,Y is closer to Left(i-1) X,Y than Right(i-1)X,Y

    for i=2:numframes;
        clear DistLL and DistLR
    %dist between  Left(i) and Left(i-1); r = sqrt((head(:,1) - targ(:,1)).^2 + (head(:,2) - targ(:,2)).^2)
    DistLL = sqrt((Left(i,1) - Left(i-1,1)).^2 + (Left(i,2) - Left(i-1,2)).^2);
    DistLR = sqrt((Left(i,1) - Right(i-1,1)).^2 + (Left(i,2) - Right(i-1,2)).^2);
   
    if DistLL < DistLR-20 & Left(i,:)==Cear1(i,:); %stay
       Left(i,:)=Cear1(i,:);
       Right(i,:)=Cear2(i,:);
    elseif DistLL < DistLR-20 & Left(i,:)==Cear2(i,:);%stay
       Left(i,:)=Cear2(i,:);
       Right(i,:)=Cear1(i,:);
    elseif DistLL > DistLR +10 & Left(i,:)==Cear1(i,:); 
       Left(i,:)=Cear2(i,:);
       Right(i,:)=Cear1(i,:); %switch
    
    elseif DistLL > DistLR +10 & Left(i,:)==Cear2(i,:); %switch
       Left(i,:)=Cear1(i,:);
       Right(i,:)=Cear2(i,:); %switch
       
       elseif DistLL > (DistLR -10) & DistLL < DistLR; %values close, click on image
        %click
        imshow(img(:,:,:,i));
        [X,Y]=getpts;
        Left(i,:)=[X(1),Y(1)];
        Right(i,:)=[X(2),Y(2)];
        
    elseif DistLL > (DistLR) & DistLL < DistLR+10; %values close, click on image
        %click
        imshow(img(:,:,:,i));
        [X,Y]=getpts;
        Left(i,:)=[X(1),Y(1)];
        Right(i,:)=[X(2),Y(2)];
    else 
        imshow(img(:,:,:,i));
        [X,Y]=getpts;
        Left(i,:)=[X(1),Y(1)];
        Right(i,:)=[X(2),Y(2)];
    end
    
    
    end
    
    
%     if DistLL < DistLR & Left(i,:)==Cear1(i,:);
%        Left(i,:)=Cear2(i,:);
%        Right(i,:)=Cear1(i,:);
%     elseif DistLL > DistLR & Left(i,:)==Cear2(i,:);
%        Left(i,:)=Cear1(i,:);
%        Right(i,:)=Cear2(i,:);
%     elseif DistLL+10 > DistLR;
%         
%         imshow(img(:,:,:,i));
%         [X,Y]=getpts;
%         Left(i,:)=[X(1),Y(1)];
%         Right(i,:)=[X(2),Y(2)];
%     end
%     figure
%     imagesc(img(:,:,:,i));hold on
%     plot(Left(1:i,1),Left(1:i,2),'g*');hold on
%     plot(Right(1:i,1),Right(1:i,2),'r*');hold on
%     
    figure
    imagesc(img(:,:,:,i));hold on
    plot(Left(1:i,1),Left(1:i,2),'g');hold on
    plot(Right(1:i,1),Right(1:i,2),'r');hold on
   %current
    plot(Left(i,1),Left(i,2),'go');hold on
    plot(Right(i,1),Right(i,2),'ro');hold on
 
% for i=1:numframes
% imshow(img(:,:,:,i));hold on
% plot(Left(i,1),Left(i,2),'c*');hold on
% plot(Right(i,1),Right(i,2),'r*');hold on
% %plot(CentroidB(i,1),CentroidB(i,2),'y*');hold on
% drawnow
% %mov(i)=getframe(gcf);
% end
    



end

        






