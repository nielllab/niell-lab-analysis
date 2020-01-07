%do for each animal, all sessions? or all animals together?
%clear all; close all
set(groot,'defaultFigureVisible','on') %disable figure plotting
savePDF=1; dbstop if error
% pname={'T:\PreyCaptureNew\Cohort3\deInterlaceTest\'};
 % pname={'T:\PreyCaptureNew\Cohort3\J463b(white)\110119\Approach\';
% pname={'T:\PreyCaptureNew\Cohort3\J463b(white)\110719\Approach\'};

deInter = 0;
doAcc = 1;

fileList=[] ;fileListR=[] ;fileListL=[] ; TSfileList=[];
%finds all files w/top.csv in the name, used for all camera files
%finds all files with acc.csv in the name for the accelerometers
for i=1:length(pname)
    fileList = [fileList; dir([pname{i} '*resnet50_Top*.csv'])];
    
    %     TSfileList = [TSfileList; dir([pname{i} '*topTS*.csv'])];
    
end


%%

for j=1% j=1:length(fileList) %%% loop over all top camera files
    if savePDF
        psfilename = 'C:\analysisPS.ps';
        if exist(psfilename,'file')==2;
            delete(psfilename);end
    else
        psfilename = 'C:\analysisPS.ps';
        
    end
    clear path Rfname Lfname %TSfileList
    %     path = fileList(j).folder
    %     TSfileList = [fileList; dir([path '*topTS*.csv'])];

    fname=fullfile(fileList(j).folder,fileList(j).name);
    
    %%% parse top camera filename
    sname = split(fname,'_');
    ani = sname{1}(end-4:end);
    sessionnum = sname{3}(end);
    date = sname{4};
    clipnum = sname{5}(~isletter(sname{5}));
    
    %%% store parameters from filename
    Data(j).ani= {ani}; %Data(j).trialtype = {trialtype};
    Data(j).sessionnum = {sessionnum};
    Data(j).date = {date};
    Data(j).clipnum = {clipnum};
    Data(j).Data = csvread(fname,3,0);
   
    %%% generate left and right eye filenames
    
    %Rfname = strrep(fname,'Top','Eye1r');
    Rfname = strrep(fname,'Top','Eye');
    if deInter
        Rfname = strrep(Rfname,'DeepCut','_DeInter2_100DeepCut');
    end
    Rfname = strrep(Rfname,'top','eye1r');
    % Rfname = strrep(Rfname,'900000','1030000');
    Rfname = strrep(Rfname,'Aug15','Jul12');
    %   Lfname = strrep(Rfname,'Eye1r','Eye2l');
    Lfname = strrep(Rfname,'eye1r','eye2l');
    
    %%% load left and right eye data
    Data(j).DataR = (csvread(Rfname,3,0))
    %         Data(j).DataR=downsamplebin(Data(j).DataR,1,2,1)
    Data(j).DataL=(csvread(Lfname,3,0));
    %         Data(j).DataL=downsamplebin(Data(j).DataL,1,2,1)
    
    %%% check length of data vs expectation
    Data(j).difTR = length(Data(j).Data)-length(Data(j).DataR)
    Data(j).difTL = length(Data(j).Data)-length(Data(j).DataL)
    Data(j).difRL = length(Data(j).DataR)-length(Data(j).DataL)
    
    
    %%% extract headpoints and cricket points from top camera
    aligned = alignHead(fname,8,0,psfilename,.90, .95)
    Data(j).mouse_xyRaw=aligned.mouse_xy;
    Data(j).mouseVRaw=aligned.mouseSp;
    Data(j).thetaRaw=aligned.theta;
    Data(j).dThetaRaw=aligned.dTheta;
    Data(j).cricketxyRaw=aligned.crick_xy;
    Data(j).cricketVRaw=aligned.crickSp;
    Data(j).cricketHRaw=aligned.crickH;
    Data(j).cricket_pHRaw=aligned.crick_pH;
    Data(j).rangeRaw=aligned.range;
    Data(j).azRaw=aligned.az;
    Data(j).cricketP=aligned.crick_p;
    Data(j).cricketThetaRaw=aligned.cricketTheta;
    Data(j).ThetaFract=aligned.ThetaFract;
    %    Data(j).dThetaFract=aligned.dThetaFract;
    %     Data(j).longTheta=aligned.longTheta;
    %     Data(j).longThetaFract=aligned.longThetaFract;
    
    %%% read in accelerometer data
    if doAcc
        accname = strcat(sname{1},'_','acc','_',sname{3},'_',date,'_',clipnum,'.dat')
        accData = getAcc(accname,psfilename);
        Data(j).accTS=accData.accTs;  %%% timestamps
        Data(j).rawAcc=accData.rawAcc;   %%% raw voltages
        Data(j).accTrace =accData.accTrace;  %%% voltages converted to m/sec2 and rad/sec
    end
    
    %%% analyze right eye points
    Data(j).xR=640 - Data(j).DataR(:,2:3:end); %%% flip right eye since this camera is reversed in bonsai
    Data(j).yR=480 - Data(j).DataR(:,3:3:end); %%% put into cartesian coords (origin lower left), instead of image coords (origin in upper left corner)
    Data(j).RLikelihood=Data(j).DataR(:,4:3:end);
    [Data(j).Rthetaraw,Data(j).Rphiraw,Data(j).EllipseParamsR,Data(j).ExtraParamsR,Data(j).goodReye, Data(j).ngoodR, Data(j).RcalR,Data(j).RcalM, Data(j).scaleR] = EyeCameraCalc1(length(Data(j).xR(:,1)), Data(j).xR,Data(j).yR, Data(j).RLikelihood,psfilename)
    Data(j).XRcentraw=Data(j).EllipseParamsR(:,1);  Data(j).YRcentraw=Data(j).EllipseParamsR(:,2);
    
    figure;imagesc(Data(j).RLikelihood); title('R eye Likelihood');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append'); end
    close(gcf)
    
    %%% analyze left eye points
    Data(j).xL= Data(j).DataL(:,2:3:end);
    Data(j).yL=480 - Data(j).DataL(:,3:3:end);    %%% put into cartesian coordinates, instead of image (origin in upper left corner)
    Data(j).LLikelihood=Data(j).DataL(:,4:3:end);
    [Data(j).Lthetaraw,Data(j).Lphiraw,Data(j).EllipseParamsL,Data(j).ExtraParamsL,Data(j).goodLeye,Data(j).ngoodL,Data(j).LcalR,Data(j).LcalM, Data(j).scaleL] = EyeCameraCalc1(length(Data(j).xL(:,1)),Data(j).xL,Data(j).yL, Data(j).LLikelihood,psfilename)
    Data(j).XLcentraw=Data(j).EllipseParamsL(:,1);  Data(j).YLcentraw=Data(j).EllipseParamsL(:,2);
    
    Data(j).RRadRaw = (Data(j).EllipseParamsR(:,3)+ Data(j).EllipseParamsR(:,4))/2;
    Data(j).LRadRaw = (Data(j).EllipseParamsL(:,3)+ Data(j).EllipseParamsL(:,4))/2;
    
    figure;imagesc(Data(j).LLikelihood); title('L eye Likelihood');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append'); end
    close(gcf)
    
    
    %fxn to adjust timing between videos of different lengths
    TSfname = strcat(ani,'_','topTS','_',sname{3},'_',date,'_',clipnum,'.csv');
    cd(fileList(j).folder)
    %read in timestamp files
    if ~isfile(TSfname)
        disp('no timestamp file');
        
    else
        
        %%% load top camera timestamps
        topTSfile=fullfile(fileList(j).folder,TSfname);
        TopTs = dlmread(topTSfile);
        TopTs= TopTs(:,1)*60*60 + TopTs(:,2)*60 + TopTs(:,3);  %%% data is read as hours, mins, secs
        startT = TopTs(1);
        Data(j).TopTs = TopTs;
        
        %%% load right eye camera timestamps
        rTSfile = strrep(topTSfile,'top','eye1r');
        RTS = dlmread(rTSfile);
        RTS= RTS(:,1)*60*60 + RTS(:,2)*60 + RTS(:,3);
        if deInter
            RTSnew = zeros(size(RTS,1)*2,1);
%             RTSnew(2:2:end) = RTS;
%             RTSnew(1:2:end) = RTS- 0.5*median(diff(RTS));
            RTSnew(2:2:end) = RTS  +0.5*median(diff(RTS)) ;
            RTSnew(1:2:end) = RTS ;
            RTS = RTSnew;
        end
        
        startR = RTS(1);
        Data(j).RTS = RTS;
        
        %%% load left eye camera timestamps
        lTSfile =  strrep(rTSfile,'eye1r','eye2l');
        LTS = dlmread(lTSfile);
        LTS= LTS(:,1)*60*60 + LTS(:,2)*60 + LTS(:,3);
        if deInter
            LTSnew = zeros(size(LTS,1)*2,1);
           % LTSnew(2:2:end) = LTS;
            % LTSnew(1:2:end) = LTS- 0.5*median(diff(LTS));
            LTSnew(2:2:end) = LTS  +0.5*median(diff(LTS)) ;
            LTSnew(1:2:end) = LTS ;
            LTS = LTSnew;
        end
        startL = LTS(1);
        Data(j).LTS = LTS;
        
        %%% determine new timestamps, based on maximum possible window
        %%% given start/stop times
        start=max([startT,startR,startL]);
        endT = min([TopTs(end-1),RTS(end),LTS(end)])  %%% TopTs goes to end-1 since dTheta only goes to end-1
        xq=(start:1/30:endT)';  %%% 30Hz resample
        
        Data(j).usedTS=xq;
        
        %%% resample top cam values
        adjustedTS = adjustTimingTop(TopTs,xq,Data(j).azRaw, Data(j).rangeRaw,Data(j).mouse_xyRaw,Data(j).mouseVRaw,...
            Data(j).cricketxyRaw,Data(j).cricketVRaw, Data(j).cricketP, Data(j).thetaRaw,Data(j).cricketHRaw,Data(j).cricket_pHRaw, Data(j).cricketThetaRaw);
        
        Data(j).az =(adjustedTS.azAdj)';
        Data(j).theta =(adjustedTS.headThetaAdj)';
        Data(j).dtheta = interp1(TopTs(1:end-1),Data(j).dThetaRaw,xq);
        Data(j).mouse_xy =adjustedTS.mousexyAdj;
        Data(j).mouseV =(adjustedTS.mouseVAdj)';
        Data(j).range=(adjustedTS.rangeAdj)';
        Data(j).cricketxy =adjustedTS.cricketxyAdj;
        Data(j).cricketV = (adjustedTS.cricketVAdj)';
        Data(j).cricketP=adjustedTS.cricketPAdj;
        Data(j).cricketHead=adjustedTS.cricketHAdj;
        Data(j).cricket_pH=adjustedTS.cricketpHAdj;
        Data(j).cricketTheta=adjustedTS.cricketThetaAdj;
        
        %%% interpolate right eye points
        Data(j).XRcent =interp1(RTS,Data(j).XRcentraw,xq);
        Data(j).YRcent =interp1(RTS,Data(j).YRcentraw,xq);
        Data(j).Rtheta =interp1(RTS,Data(j).Rthetaraw,xq);
        Data(j).Rphi = interp1(RTS,Data(j).Rphiraw,xq);
        Data(j).RRad = interp1(RTS,Data(j).RRadRaw,xq);
       
        %%% interpolate left eye points
        Data(j).XLcent =interp1(LTS,Data(j).XLcentraw,xq);
        Data(j).YLcent =interp1(LTS,Data(j).YLcentraw,xq);
        Data(j).Ltheta =interp1(LTS,Data(j).Lthetaraw,xq);
        Data(j).Lphi =interp1(LTS,Data(j).Lphiraw,xq);
        Data(j).LRad = interp1(LTS,Data(j).LRadRaw,xq);
        
        figure;subplot(1,2,1)
        plot(Data(j).Rtheta); title('R eye Theta');
        subplot(1,2,2)
        plot(Data(j).Ltheta); title('L eye Theta');
        
        if savePDF,set(gcf, 'PaperPositionMode', 'auto');print('-dpsc',psfilename,'-append');   end
        close(gcf)
        
        figure;plot(Data(j).RRad); hold on; plot(Data(j).LRad);
        title('interp pupil radius');
        
        %%% interpolate accelerometers
        if doAcc
            
            %%% 
            for i = 1:6
                Data(j).accResamp(:,i) = interp1(Data(j).accTS,Data(j).accTrace(:,i),xq);
                Data(j).rawAccResamp(:,i)= interp1(Data(j).accTS,Data(j).rawAcc(:,i),xq)
            end
            
            %%% there is a constant offset in the accelerometer timestamps coming in through labjack.
            %%% we calculate this by taking xcorr of head rotation from DLC 
            %%% and head rotation from gyros (should be identical) and
            %%% we then shift the data accordingly.
            gyro3=(Data(j).accResamp(:,6)-nanmean(Data(j).accResamp(:,6))); % gyro 3 = yaw
            dth=Data(j).dtheta;
            [corr lags]=nanxcorr(gyro3,dth,100,'coeff');
            
            figure;  plot(lags,corr); axis square; title('acc yaw, d head th')
            [mx ind] = max(corr);
            drift = lags(ind);
            Data(j).accShift = circshift(Data(j).accResamp,-drift,1);
            Data(j).accXcorrMax = mx;
            Data(j).accXcorrLag = drift;
            Data(j).rawAccShift = circshift(Data(j).rawAccResamp,-drift,1)          
        end
        
        %%% calculate differences
        Data(j).dxR = diff(Data(j).XRcent);
        Data(j).dxL = diff(Data(j).XLcent);
        Data(j).dyR = diff(Data(j).YRcent);
        Data(j).dyL = diff(Data(j).YLcent);
        Data(j).dxRTheta = diff(Data(j).Rtheta);
        Data(j).dxLTheta = diff(Data(j).Ltheta);
        Data(j).dxRPhi = diff(Data(j).Rphi);
        Data(j).dxLPhi = diff(Data(j).Lphi);
        Data(j).dth = Data(j).dtheta;


        %%% a few figures ...
           
        figure
        plot(xcorr(Data(j).dxRTheta,Data(j).dth(1:end-1),30,'coeff'))
        hold on
        plot(xcorr(Data(j).dxLTheta,Data(j).dth(1:end-1),30,'coeff'))
        plot(xcorr(Data(j).dxRTheta,-Data(j).dxLTheta,30,'coeff'))
        legend('R','L','R-L')
        title('eye theta & head theta');
        %         if savePDF
        %             set(gcf, 'PaperPositionMode', 'auto');
        %             print('-dpsc',psfilename,'-append');
        %         end
        %         close(gcf)
        
        figure
        plot(xcorr(Data(j).dxRPhi,Data(j).dth(1:end-1),30,'coeff'))
        hold on
        plot(xcorr(Data(j).dxLPhi,Data(j).dth(1:end-1),30,'coeff'))
        plot(xcorr(Data(j).dxRPhi,Data(j).dxLPhi,30,'coeff'))
        legend('R','L','R-L')
        title('eye phi & head theta');
        %         if savePDF
        %             set(gcf, 'PaperPositionMode', 'auto');
        %             print('-dpsc',psfilename,'-append');
        %         end
        %         close(gcf)
        %
        
        
    end
    if savePDF
        pSname='T:\PreyCaptureAnalysis\Data\singleVid_pdfs\';
        C={ani, date, sessionnum, clipnum};
        filen=sprintf('%s%s%s%s',ani,date,sessionnum,clipnum,'originalAcc','.pdf')
        pdfilename=fullfile(pSname,filen)
        dos(['ps2pdf ' psfilename ' ' pdfilename]);
        delete(psfilename);
    else
        pFile='T:\PreyCaptureAnalysis\Data\';
        
    end
    
    
end
pFile='T:\PreyCaptureAnalysis\Data\';
 afilename=sprintf('%s',ani,'_originalACC','.mat')
%afilename=sprintf('%s',ani,'_121019','.mat')

save(fullfile(pFile, afilename))
%save('J463c_test_data.mat')

