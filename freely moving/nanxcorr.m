function [cc lags] = nanxcorr(x,y,maxlag,normalization);
%%% calculates xcorr ignoring NaNs without altering timing
%%%
%%% achieves this by manually shifting y, then removing missing data and
%%% doing auto-correlation, which doesn't depend on relative timing
%%%
%%% uses same inputs/outputs as standard xcorr
%%% except normalization option = 'zero' does mean subtraction before coeff normalization
%%%
%%% cmn 08/2019

%%% check inputs
if ~exist('normalization','var')
    normalization = 'coeff';
end
if ~exist('maxlag','var')
    maxlag = 25;
end

if strcmp(normalization, 'zero')
    x = x-nanmean(x);
    y = y-nanmean(y);
    normalization = 'coeff';
end


%%% inputs must be vectors, of same orientation (so we make them rows)
if ~isvector(x) | ~isvector(y)
    display('must be vectors');
    return
end
if ~isrow(x), x = x'; end
if ~isrow(y), y = y'; end

%%% calculate correlation across lags
lags = -maxlag:maxlag;
for i = 1:length(lags);
    yshift = circshift(y,lags(i));
    use = ~isnan(x+yshift);
    cc(i,1) = xcorr(x(use),yshift(use),0,normalization);
end
