function f =plotFiltSections(sections,anatomy,x,y,histSection,type, range,cmap,mapRange,sigma)

if ~exist('range','var') | isempty(range)
    range = 1:length(sections)
    use_subplot=1;
else
    use_subplot=0;
end

f=figure;

size(type)
size(histSection)


x=x+500;
y=y+500;



for i = range

   
    sectionx = sections(i).coords(:,1)+500;
     sectiony = sections(i).coords(:,2)+500;
    hold on
    if use_subplot
        subplot(2,3,i)
    end
    
    %sigma =100;
    f= fspecial('gaussian',500,sigma);
    f = f/max(f(:));
    i
    size(histSection)
    cells = find(histSection==i);
    
    space = zeros(1000,1000);
    sites = histSection==i & ~isnan(type);
    space(sub2ind(size(space),x(sites),y(sites)))=type(sites);
    truefilter = imfilter(space,f,'same');   
%     figure
%     imagesc(truefilter')
%     axis xy
    
    allspace = zeros(1000,1000);
    sites = histSection==i & ~isnan(type);
    space(sub2ind(size(space),x(sites),y(sites)))=1;
    allfilter = imfilter(space,f,'same');
    allfilter(allfilter<0.1)=0;
%     figure
%     imagesc(allfilter')
%     axis xy
    

    inside = zeros(1000,1000);
    for xind = min(sectionx):max(sectionx);
        sectiony(sectionx==xind);
        inside(xind,min(sectiony(sectionx==xind)):max(sectiony(sectionx==xind)))=1;
    end
 
    
    prob = truefilter./allfilter;
    prob(prob==Inf)=nan;
    prob(~inside)=nan;
 
    h=imagesc(prob',mapRange);
    set(h,'alphadata',~isnan(prob'))
    axis xy
    colormap(cmap);
    hold on
    plot(sectionx,sectiony ,'b.','MarkerSize',2)
    axis([0 1000 0 1000])
    axis square
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
     
end
colordef white