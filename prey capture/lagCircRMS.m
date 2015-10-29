function xc = lagCircRMS(x,y,lags);
for lag = 1:length(lags)
   % xc(lag) = mean((mod(circshift(y,lags(lag))-x+pi,2*pi)-pi).^2);
   xc(lag) = mean(cos(x-circshift(y,lags(lag))));
end
