function cv = calcCVosi(tuning);
if min(tuning)<0
    tuning = tuning-min(tuning);
end
cv = sum(tuning.*exp((1:length(tuning))*sqrt(-1)*2*pi/(0.5*length(tuning))));
cv = cv/sum(tuning);
