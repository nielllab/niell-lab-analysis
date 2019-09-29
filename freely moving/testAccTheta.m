
load('accTest_J462a_090519_1_3.mat');
load('J462aAnalyzed_oneSess_090519.mat');

i = 2

%%% need to recalculate the timepoints (xq) that were used for interpolation
%%% better is probably to save this out for each session in loadAllCsv

start=max([Data(i).TopTs(1),Data(i).RTS(1),Data(i).LTS(1)]);
endT = min([Data(i).TopTs(end-1),Data(i).RTS(end),Data(i).LTS(end)])  %%% TopTs goes to end-1 since dTheta only goes to end-1
xq=(start:1/30:endT)';

%%% resample accelerometer data
accResamp = interp1(accTs,accTrace,xq);

%%% get head dtheta, and clean up large artifact values
dth = dTheta{i};
dth(abs(dth)>30)=NaN;

drift = 262; %%% I had to identify this by eye ... increasing value until vergence and acc2 lined up


for j = 4:6; %%%%loop over acceleratometers (6 = yaw)    
    figure
    plot(dth);
    hold on
    plot((circshift(accResamp(:,j),drift) - nanmean(accResamp(:,j)))*3);
    % ylim([-40 40])
    %             plot(xcorr(accResamp(good,j),dth(good),300,'coeff'))
    
    legend('dtheta','acc ch6')
end


%%% calculate vergence
verg = thetaL{2} - thetaR{2};

%%% loop over all accelerometer channels
for j = 1:3
    figure
    subplot(2,1,1);
    plot(thetaL{2} - thetaR{2}); ylim([-90 90])
    subplot(2,1,2);
    plot(circshift(tilt,drift))
end

%%% channel 2 is tilt, do some smoothing since high freq is also linear acceleration
tilt = conv(accResamp(:,2),[1 1 1 1 1],'same');
figure
plot(verg); hold on;
plot(circshift(tilt*4-20,262)); ylim([-20  90])
legend('vergence', 'acc ch2 (tilt)')

%%% channel 1 appears to be roll, matches with left-right phi difference.
roll = conv(accResamp(:,1),[1 1 1 1 1],'same');
eyeRot = phiR{2} - phiL{2};
figure
plot(eyeRot); hold on;
plot(circshift(roll*6-60,262)); ylim([-60  60])
legend('phi diff', 'acc ch1 (roll)')

