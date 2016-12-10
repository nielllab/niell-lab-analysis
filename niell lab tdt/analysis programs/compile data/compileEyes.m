clear all
close all 

dbstop if error

batchEyetracking_angie; %%% load batch file
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
use =  find(strcmp({files.notes},'good data'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;

for i = 1:length(use)
    
    cfile = [{[pathname '\' files(use(i)).dir '\' files(use(i)).prehigh_camera '.mat']}; {[pathname '\' files(use(i)).dir '\' files(use(i)).posthigh_camera '.mat']};
        {[pathname '\' files(use(i)).dir '\' files(use(i)).prelow_camera '.mat']};{[pathname '\' files(use(i)).dir '\' files(use(i)).postlow_camera '.mat']}]';

%   % sessionNum(cellrange)=i;
%     for j=1:length(cellrange); expt{cellrange(j)} = files(use(i)).expt; end
%     if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline, end;
%     if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi, end;
%     if strcmp(files(use(i)).treatment,'5HT'), treatment(cellrange)=ht, end;
%     if strcmp(files(use(i)).treatment,'Ketanserin'), treatment(cellrange)=ketanserin, end;
%     if strcmp(files(use(i)).treatment,'KetanserinDOI'), treatment(cellrange)=ketandoi, end;
%     if strcmp(files(use(i)).treatment,'MGluR2'), treatment(cellrange)=mglur2, end;
%     if strcmp(files(use(i)).treatment,'MGluR2DOI'), treatment(cellrange)=mglur2doi, end;
%     if strcmp(files(use(i)).treatment,'Lisuride'), treatment(cellrange)=lisuride, end;


%plotting fine, but saving values as NANs???
if ~isempty(files(use(i)).blockHigh{1}) & ~isempty(files(use(i)).blockHigh{2})
    
    for prepost =1:2
        eyes = getEyes_angie(cfile{:,prepost}, files(use(i)).blockHigh{prepost},1);
        radHigh{:,prepost} = eyes.rad
        rInterpHigh{:,prepost}=eyes.rInterp
        vInterpHigh{:,prepost}=eyes.vInterp
        tHigh{:,prepost} = eyes.t
    end
end   
   if ~isempty(files(use(i)).blockLow{1}) & ~isempty(files(use(i)).blockLow{2})

   for prepost =1:2
        eyes = getEyes_angie(cfile{:,prepost}, files(use(i)).blockHigh{prepost},1);
        radLow{:,prepost} = eyes.rad
        rInterpLow{:,prepost}=eyes.rInterp
        vInterpLow{:,prepost}=eyes.vInterp
        tLow{:,prepost} = eyes.t
   end
end
end


figure
subplot(2,2,1)
set(gcf,'Name','high contrast prepost');
plot(rInterpHigh{1});hold on; plot(vInterpHigh{1},'g'); ylim([0 30]); axis xy;title('Pre')
subplot(2,2,2)
plot(rInterpHigh{2}); hold on; plot(vInterpHigh{2},'g');axis xy;title('Post');
legend('radius','velocity')
subplot(2,2,3)
plot(rInterpHigh{1},vInterpHigh{1},'.'); xlim([10 30]);ylim([0 40]); 
hold on;lsline;ylabel('velocity cm/sec');xlabel('radius'); title('Pre')
subplot(2,2,4)
plot(rInterpHigh{2},vInterpHigh{2},'.'); xlim([10 30]);ylim([0 40]);
hold on;lsline;ylabel('velocity cm/sec');xlabel('radius'); title('Post')

figure
subplot(2,2,1)
set(gcf,'Name','low contrast prepost');
plot(rInterpLow{1});hold on; plot(vInterpLow{1},'g');  ylim([0 30]);axis xy
title('Pre');
subplot(2,2,2)
plot(rInterpLow{2}); hold on; plot(vInterpLow{2},'g'); ylim([0 30]);axis xy
legend('radius','velocity'); title('Post');
subplot(2,2,3)
plot(rInterpLow{1},vInterpLow{1},'.'); xlim([10 30]); ylim([0 40]);hold on; axis xy;
lsline;ylabel('velocity cm/sec');xlabel('radius'); title('Pre')
subplot(2,2,4)
plot(rInterpLow{2},vInterpLow{2},'.');xlim([10 30]); ylim([0 30]);hold on; axis xy;
lsline;ylabel('velocity cm/sec');xlabel('radius');title('Post')


