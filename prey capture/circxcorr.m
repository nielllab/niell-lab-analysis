function xc = circxcorr(x,y,lags);
for lag = 1:length(lags)
    xc(lag) = circ_corrcc(x,circshift(y,lags(lag)));
end