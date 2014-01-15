function f =plotSections(sections,anatomy,histox,histoy,histSection,labels, range)

colordef black
if ~exist('range','var')
    range = 1:length(sections)
    use_subplot=1;
else
      use_subplot=0;
end

f=figure;

for i = range
    hold on
    if use_subplot
        subplot(2,3,i)
    end
    plot(sections(i).coords(:,1),sections(i).coords(:,2) ,'b.','MarkerSize',2)
    axis([-500 500 -500 500])
    axis square
    hold on
%     for j=1:length(anatomy)
%         if anatomy(j).section==i
%             plot(anatomy(j).siteXY(1,:),anatomy(j).siteXY(2,:),'g.','MarkerSize',10);
%         end
%     end
    %axis off
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    cells = find(histSection==i);
    if ~isempty(cells)
        scatter(histox(cells)+30*(rand(size(cells))-0.5),histoy(cells),[],labels(cells,:),'.');
    end
end
colordef white