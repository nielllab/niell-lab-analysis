function y = interpNan(x,r,method);
%%% interpolates across NaNs, but only for runs less than length r
%%% leaves NaNs at start/end since we don't want to extrapolate
%%% takes option input 'method' to specify method to pass to interp1
if ~exist('method','var')
    method = 'linear';
end

%%% perform initial interpolation
y = interp1(find(~isnan(x)),x(~isnan(x)),1:length(x),method);

%%% find start and end of NaN runs
bad = isnan(x);
starts = find(diff(bad)>0) + 1;

if isnan(x(1))
    try             %starts can be either row or column, depending on input
        starts = [1 starts];
    catch
        starts = [1 starts'];
    end
end

ends = find(diff(bad)<0);

%%% leave runs > r as bad, mark runs < r as not bad
for i = 1:length(ends)  %%% note, there can be 1 more starts than ends, if the sequence ends with a NaN. but we want to leave that in
    if (ends(i)-starts(i))<r & starts(i)~=1   %%% don't fill in NaNs at beginning
        bad(starts(i):ends(i))=0;
    end
end

%%% put back in NaNs for bad interpolations
y(bad)=NaN;

