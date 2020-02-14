%%% test program to determine best way for lining up accelerometer data
%%% given offset in timestamps (resulting from labjack?)
%%% old method did standard interpolation using incorrect points, 
%%% but then xcorr to determine a shift (ergo must be multiple of framerate)
%%% new way tests a range of offset to add for the timestamp
%%% both compare interpolated acc data to head dtheta data

for i = 1:length(Data);
    accPre = Data(i).accTrace(:,6);
    ts = Data(i).accTS;
    xq = Data(i).usedTS;
    
    dth=diff(Data(i).theta); dth(dth<-pi)= dth(dth<-pi)+2*pi; dth(dth>pi) = dth(dth>pi) -2*pi;
    
    offsets = -10:0.005:10;
    for j = 1:length(offsets)
        newinterp = interp1(ts+offsets(j),accPre, xq);
        xc(j)= nanxcorr(newinterp(1:end-1),dth,0, 'coeff');
    end
    
    [max_xcorr max_ind] = max(xc);
    accPost = interp1(ts+offsets(max_ind),accPre, xq);
    
    
    eyes = 0.5*(Data(i).dxRTheta + Data(i).dxLTheta);
    
    figure
    subplot(1,2,1)
    
    plot(Data(i).accShift(1:end-1,6),eyes,'.');
    axis equal; axis([-20 20 -20 20])
    title('shifted only')
    
    subplot(1,2,2)
    plot(accPost(1:end-1),eyes,'.');
    axis equal; axis([-20 20 -20 20])
    title('interpolated')
    
end
