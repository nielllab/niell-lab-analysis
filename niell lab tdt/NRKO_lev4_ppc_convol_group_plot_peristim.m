function [] = NRKO_lev3_ppc_convol_group_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

% take the globals from the info script that has been ran at the start
global info
global outputDir
output  = 'lev3_ppc_convol_group_peristim';
load(fullfile(outputDir,output,output), 'group')

% do a random effects analysis;
% first delete all cells with <50 spikes, used this in most papers
iMethod = 1;
keep = group.dof{iMethod}>50;
keep = keep(:,end);
iMethod = 1;
ppc = group.ppc{iMethod}(keep,:);
genotype = group.genotype(keep);
animalid = group.animalid(keep);
animals  = unique(animalid);
meanAnimal = []; varAnimal = []; genotypeAnimal = [];
for iAnimal = 1:length(animals)
  animalSel = animalid==animals(iAnimal);
  meanAnimal(iAnimal,:)   = nanmean(ppc(animalSel,:));
  varAnimal(iAnimal,:)    = nansem(ppc(animalSel,:)).^2;
  genotypeAnimal(iAnimal) = unique(genotype(animalSel));
end

% get the means according to random effects analysis
gCnt = 0; mnGG = []; smGG = [];
for iG = 0:1
  gCnt = gCnt + 1;
  sl = genotypeAnimal==iG;
  w  = 1./varAnimal(sl,:);
  mnG = repmat(sum(meanAnimal(sl,:).*w)./sum(w),[sum(sl) 1]);
  Q  = sum(w.*(meanAnimal(sl,:)-mnG).^2) - (sum(sl)-1);
  nrm = sum(w) - sum(w.^2)./sum(w);
  tausq = Q./nrm;
  tausq(tausq<0) = 0;
  tausq = repmat(tausq, [sum(sl) 1]);
  w = 1./(varAnimal(sl,:) + tausq);
  mnGG(gCnt,:) = sum(w.*meanAnimal(sl,:))./sum(w);
  smGG(gCnt,:) = sqrt(1./sum(w));
end
  
figure, errorbar(1:3,mnGG(1,:),smGG(1,:),smGG(1,:)), hold on, errorbar(1:3,mnGG(2,:),smGG(2,:),smGG(2,:),'r')





