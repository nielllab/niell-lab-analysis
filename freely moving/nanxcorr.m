function [cc lags] = nanxcorr(x,y,maxlag,normalization);
%%% calculates xcorr ignoring NaNs without altering timing
%%%
%%% achieves this by manually shifting y, then removing missing data and
%%% doing auto-correlation, which doesn't depend on relative timing
%%%
%%% uses same inputs/outputs as standard xcorr
%%%
%%% cmn 08/2019


if ~exist('normalization','var')
    normalization = 'coeff';
end

if ~exist('maxlag','var')
    maxlag = 25;
end

lags = -maxlag:maxlag;
for i = 1:length(lags);
    yshift = circshift(y,lags(i));
    use = ~isnan(x+yshift);
    cc(i,1) = xcorr(x(use),yshift(use),0,normalization);
end
