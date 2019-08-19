function adjustedData = adjustTiming(EyeTS,desiredTimestamps,theta,phi,xcent,ycent,psfilename)
if exist('psfilename','var')
    savePDF=1;
end

%adj = interp1(rawData,csvEyeTimestamps,desiredTimestamps, vars to adjust)
thetaAdj = interp1(EyeTS,theta,desiredTimestamps);
phiAdj = interp1(EyeTS,phi,desiredTimestamps);
x_centAdj = interp1(EyeTS,xcent,desiredTimestamps)
y_centAdj = interp1(EyeTS,ycent,desiredTimestamps)


figure;
subplot(2,2,1)
plot(theta,'b'); hold on; plot(thetaAdj,'g');axis square; title('theta')
subplot(2,2,2)
plot(phi,'b'); hold on; plot(phiAdj,'g'); axis square;title('phi')
subplot(2,2,3)
plot(xcent,'b');hold on; plot(x_centAdj,'g');axis square;title('x cent')
subplot(2,2,4)
plot(ycent,'b');hold on; plot(y_centAdj,'g');axis square;title('y cent')

 if savePDF
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpsc',psfilename,'-append');
    end
   close(gcf)

adjustedData.thetaAdj=thetaAdj;
adjustedData.phiAdj=phiAdj;
adjustedData.x_centAdj=x_centAdj;
adjustedData.y_centAdj=y_centAdj;

end
