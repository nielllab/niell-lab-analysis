function [] = NRKO_lev3_ppc_convol_group_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev2_ppc_convol_peristim';
output  = 'lev3_ppc_convol_group_peristim';

% loop over the various files
nDirs = length(info);
methods  = {'ppc0', 'ppc1', 'ppc2', 'ang', 'ral'};
nMethods = length(methods); 
[dofAllCat, ppcAllCat] = deal(cell(1,nMethods));
genotype = []; animalid = [];
for iDir = 1:nDirs
  
  load(fullfile(outputDir, input, info(iDir).dataname, input));  
  for iMethod = 1:nMethods
    ppcCat = cat(1,stat(iMethod).static(:).(methods{iMethod}));
    dofCat = cat(1,stat(iMethod).static(:).nspikes);    
    ppcAllCat{iMethod} = [ppcAllCat{iMethod};ppcCat];
    dofAllCat{iMethod} = [dofAllCat{iMethod};dofCat];
    nCells   = size(ppcCat,1);
    if iMethod==1
      genotype = [genotype; ones(nCells,1)*info(iDir).genotype];
      animalid = [animalid; ones(nCells,1)*iDir];
    end
  end
end


% gather the results
group.genotype = genotype;
group.animalid = animalid;
group.ppc      = ppcAllCat;
group.dof      = dofAllCat;
group.methods  = methods;
group.genotypelabels = {'0=control, 1=knockout'};

% save the data to the disk
%keyboard
mkdir(fullfile(outputDir, output));
save(fullfile(outputDir,output,output), 'group')





