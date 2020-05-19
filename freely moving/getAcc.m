function data = getAcc(accname,framerate,psfilename)
filterWin = 5; %%% timepoints to median filter over for acceleration

data=dlmread(accname)

accTs = data(:,1); %timestamps
% dname = split(accname,'_');
% expdate=dname{4};
% if sum(expdate>110319)<6
%     accTs = mod(accTs-7*60*60,24*60*60); 
%     display('pre dls')
% else
accTs = mod(accTs-8*60*60,24*60*60); %%% time is elapsed secs since midnight 1904 GMT; subtract 8 hrs to get local time (but what about daylight savings change!)
% end

figure
plot(diff(accTs))
ylabel('secs'); title(sprintf('diff of acc timestamps; median = %0.3f',median(diff(accTs))));
if exist('psfilename','var')
    savePDF=1;
end

rawAcc=data(:,2:7);
filtAcc = (medfilt1(rawAcc(:,1:3),filterWin)-2.5)*1.6; %  acc ch conversion and filtering
filtAcc(filtAcc>1)=1; filtAcc(filtAcc<-1)=-1;
acc = asind(filtAcc);

% good = abs(acc)<1;
% 
% acc(good) = real(acc(good));
% acc(~good) = NaN;
acc(:,2)=-acc(:,2); %flip sign of channel 2 & 3
acc(:,3)=-acc(:,3); %flip sign of channel 2 & 3
gyro = (data(:,5:7)-2.5)*(400/framerate); 
%convert voltage from gyro channels to degrees: Full range (5V) is set to  +/- 1000deg/sec. So ...
%5V = 2000deg/sec, 1V = 400deg/sec, so 1V = (400/framerate) deg/frame.
%So multiply gyro data by 400/framerate

%%% get data (columns 2-7)

figure
plot(accTs-accTs(1),acc);
plot(accTs-accTs(1),gyro);
xlabel('secs'); ylim([0 5]);

if exist('psfilename','var')
    savePDF=1;
    
clear data
data.accTrace(:,1:3) = acc;
data.accTrace(:,4:6) = gyro;
data.accTs = accTs;
data.rawAcc=rawAcc;

end