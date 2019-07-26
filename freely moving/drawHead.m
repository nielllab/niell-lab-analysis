function drawHead(pts,sm)
if ~exist('sm','var')  % option for 'small' figure, with no circles
    sm=0;
end
plot(pts([3 2 5 6],1),pts([3 2 5 6],2),'LineWidth',2);
hold on
plot(pts([4 8 7],1), pts([4 8 7],2),'LineWidth',2);
plot(pts([1 8],1),pts([1 8],2),'g','LineWidth',2)

if ~sm %%% 
plot(pts(1,1),pts(1,2),'g*','Markersize',12)
plot(pts([3 4 6 7],1), pts([3 4 6 7],2),'ko','Markersize',12)
plot(mean(pts(:,1)),mean(pts(:,2)),'x','Markersize',12)
end