function adjustedData = adjustTimingTop(TopTS,desiredTimestamps,az,range,mousexy,mouseV,cricketxy,cricketV,cricketP,headtheta)

azAdj = interp1(TopTS,az,desiredTimestamps);
rangeAdj = interp1(TopTS,range,desiredTimestamps);
mouseVAdj = interp1(TopTS(1:end-1),mouseV,desiredTimestamps);
cricketVAdj = interp1(TopTS(1:end-1), cricketV, desiredTimestamps);
cricketPAdj = interp1(TopTS, cricketP, desiredTimestamps);
headThetaAdj = interp1(TopTS, headtheta, desiredTimestamps);

for xy=1:2
mousexyAdj (xy,:) = interp1(TopTS,squeeze(mousexy(xy,:)),desiredTimestamps);
cricketxyAdj(xy,:) = interp1(TopTS,squeeze(cricketxy(xy,:)),desiredTimestamps);
end

adjustedData.rangeAdj=rangeAdj;
adjustedData.mouseVAdj=mouseVAdj;
adjustedData.azAdj=azAdj;
adjustedData.cricketPAdj=cricketPAdj;
adjustedData.cricketVAdj=cricketVAdj;
adjustedData.headThetaAdj=headThetaAdj;
adjustedData.mousexyAdj=mousexyAdj;
adjustedData.cricketxyAdj=cricketxyAdj;




end