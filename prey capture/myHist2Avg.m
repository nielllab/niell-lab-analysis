function [h,avgY] = myHist2Avg(x,y,xbins,ybins);
%%% 2D histogram
%%% specifies bins by centers, as does hist
%%% unlike hist2, which works like histc and has unpredictable behavior
%%% (for example, hist2 failes if ybins > xbins

for i = 1:length(xbins);
    if i ==1
        data = find(x<=mean(xbins(1:2)));
    elseif i == length(xbins)        
        data = find(x>mean(xbins(end-1:end)));
    else
        data = find(x>mean(xbins(i-1:i)) & x<=mean(xbins(i:i+1)));
    end
    h(:,i) = hist(y(data),ybins);
    avgY(i)=nanmean(y(data));
%     std(:,i)=nanstd(y(data));
%     SEm(:,i)=std(:,i)/sqrt(length(y(data)));
end
