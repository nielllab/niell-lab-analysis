function cv = calcCVosi(tuning);
if min(tuning)<0
    tuning = tuning-min(tuning);
end
cv = sum(tuning.*exp((1:8)*sqrt(-1)*2*pi/4));
cv = cv/sum(tuning);
