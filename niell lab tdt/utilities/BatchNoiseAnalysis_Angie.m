dbstop if error
clear all
close all 

batchDOIephys_filtered;
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%use =  find(strcmp({files.notes},'good data')& ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) )
%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )
use =  find(strcmp({files.notes},'good data')  & ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) & strcmp({files.expt},'102516'))
parpool
 for i = 1:length(use)
     close all
     afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
     clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
     
     movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
     load(movieFile);
   if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
         
         for prepost = 1:2
             %  noise_analysis_angie(clustfile,afile,movieFile,block,movietype, redo,periofFrames,stim_eye)
             noise_analysis_angie(clustfile,afile,movieFile,files(use(i)).blockWn{prepost},1,0,300,1)
             fit2dgabor_angie(afile, files(use(i)).blockWn{prepost},1)
         end
     end
% movietype==fl_noise==2; 
% movietype==cm_noise==1 
% movietype==mv_noise==3
%    toc 
 end
   delete(gcp('nocreate'))