clear all
close all

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})  & strcmp({files.expt},'120115') )

use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

sprintf('%d selected sessions',length(use))

prespikes = []; %%% initialize data to be saved
postspikes = [];

for i = 1:length(use)
    afile = files(use(i)).analysisfile;
    clustfile = files(use(i)).clusterfile;
    
    block = files(use(i)).predark;   %%% run for pre
    darknessAnalysis
    %%% save out any important data
    prespikes{end+1:end+length(blockSpike)} = blockSpike;
    
    
    block = files(use(i)).postdark;  %%% run for post
    darknessAnalysis
    %%% save out any important data
    postspikes{end+1:end+length(blockSpike)} = blockSpike;
end
