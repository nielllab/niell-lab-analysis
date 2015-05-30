function [xc tbase]= spikeLFPcoh(sp,t,lfp,dt,range);
tbins =0:dt:t(end);
sp_hist = hist(sp,tbins);
lfp_int = interp1(t,lfp,tbins);
shifts = round((-range/dt):(range/dt));
tbase = shifts*dt;
for i = 1:length(shifts);
   
xc(i) = sp_hist*circshift(lfp_int',shifts(i));
end
xc = xc/length(sp);

