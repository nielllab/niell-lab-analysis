function [] = NRKO_lev3_ppc_convol_group_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

%% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev2_tuning_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev3_tuningstats_group_peristim';

% loop over the various files
nDirs = length(info);
latencies = {'prestim', 'stim'};
statelabels{1} = 'move';
statelabels{2} = 'sit';
statelabels{3} = 'early_sit'; % <10 s.
statelabels{4} = 'late sit';
nStates  = length(statelabels);
nLatencies = length(latencies);
genotype = []; animalid = [];
dimord = 'state_latency';
%%
totUnits = 0;
clear ppcAllCat
for iDir = 1:nDirs
    iDir
    clear dat data
    load(fullfile(outputDir, input, info(iDir).dataname, input));  
    inp2 = load(info(iDir).datafile);
  
    nUnits = length(inp2.unitdataSession)
    nUnits = length(Stat.osi)
    for iUnit = 1:nUnits
        unitCnt =  totUnits + iUnit;                   
        unitinfo(unitCnt,1) = inp2.unitdataSession{iUnit}.GT*0+info(iDir).genotype;
        unitinfo(unitCnt,2) = inp2.unitdataSession{iUnit}.pinp;
        unitinfo(unitCnt,3) = inp2.unitdataSession{iUnit}.responsive;
        unitinfo(unitCnt,4) = inp2.unitdataSession{iUnit}.layer;
        unitinfo(unitCnt,5) = inp2.unitdataSession{iUnit}.inhibitory;             
        if strcmp(info(iDir).stimType,'drift')
            unitinfo(unitCnt,6) = 1;  
        else
            unitinfo(unitCnt,6) = 2;
        end

        osi_cat(unitCnt) = Stat.osi(iUnit);
        loco_cat(unitCnt,:) = Stat.modIndx_locomotion(:,iUnit);
        if length(Stat.modIndx_stim)==nUnits
            stim_cat(unitCnt) = Stat.modIndx_stim(iUnit);            
        else
            stim_cat(unitCnt) = NaN;
        end
        rate_cat(unitCnt,:,:) = Stat.mnRate(:,:,iUnit);
        animalinfo{unitCnt} = inp2.data.tank;            
    end
    totUnits = totUnits + nUnits;   
end
%%
rate_stat.unitinfo = unitinfo;
rate_stat.animalinfo = animalinfo;
rate_stat.osi = osi_cat;
rate_stat.stim_mod = stim_cat;
rate_stat.loco_mod = loco_cat;
rate_stat.rate = rate_cat;
rate_stat.unitinfo_dim = {'wt_N2B_N2A', 'pinp_nopinp', 'responsive', 'layer', 'inhibitory', 'drift1_bar2'};
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
save(filename,'rate_stat');
