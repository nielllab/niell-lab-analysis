clear all
close all

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})  & strcmp({files.expt},'120115') )

use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

sprintf('%d selected sessions',length(use))

for i = 1:length(use)
    afile = files(use(i)).analysisfile;
    clustfile = files(use(i)).clusterfile;
    
    block = files(use(i)).prewn;   %%% run for pre
    analyzeDarkness;
    %%% save out any important data
    
    block = files(use(i)).postwn;  %%% run for post
    analyzeDarkness;
    %%% save out any important data
end
