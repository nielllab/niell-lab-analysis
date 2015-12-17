function [stat] = Lur_statistics(cfg, data, individuals, conditions, conditionlabels)

if ~isfield(cfg, 'variance_estimator'), cfg.variance_estimator = 'pooled'; end
if ~isfield(cfg, 'tausq_est'), cfg.tausq_est = 'weighted'; end % else unweighted
if ~isfield(cfg, 'semiweighted_adaptive'), cfg.semiweighted_adaptive = 'no'; end % else unweighted
if ~isfield(cfg, 'semiweighted_adaptive_niter'), cfg.semiweighted_adaptive_niter = 1000; end
% options
% 
% test data
% four animals
% conditions = [ones(1,100)*1 ones(1,100)*2];
% data1 = [];
% individuals = [];
% for k = 1:10
%   data1 = [data1 randn(1,10)+2*rand];
%   individuals  = [individuals ones(1,10)*k];
% end
% data2= [];
% for k = 1:10
%   data2 = [data2 1+randn(1,10)+2*rand];
%   individuals  = [individuals ones(1,10)*(k+10)];
% end
% data = [data1 data2];
% 
% remove nans from the data and create row vectors
todel       = ~isnan(data);
data        = data(todel);
data        = data(:);
conditions  = conditions(todel);
conditions  = conditions(:);
individuals = individuals(todel);
individuals = individuals(:);

% get the unique conditions
uniqueConditions = unique(conditions);
nConditions      = length(uniqueConditions);

%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE POOLED OR WEIGHTD ESTIMATOR
[mnCondition_weighted,smCondition_weighted] = deal([]);
for iCondition = 1:nConditions
  isCondition = conditions==uniqueConditions(iCondition);
  mnCondition_weighted(iCondition) = mean(data(isCondition));
  smCondition_weighted(iCondition) = sem(data(isCondition));
end

%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE THE UNWEIGHTED MEANS ACROSS INDIVIDUALS
;
[mnCondition_unweighted, smCondition_unweighted] = deal([]);
for iCondition = 1:nConditions
  
  hasCondition = conditions==uniqueConditions(iCondition);  
  uniqueIndividuals = unique(individuals(hasCondition));
  nIndividuals      = length(uniqueIndividuals);
  [mnIndividual,smIndividual, varIndividual,dofIndividual,conditionIndividual] = deal([]);
  for iIndividual   = 1:nIndividuals
    isIndividual = individuals==(uniqueIndividuals(iIndividual)); 
    mnIndividual(iIndividual) = mean(data(isIndividual&hasCondition));
    smIndividual(iIndividual) = sem(data(isIndividual&hasCondition));
    varIndividual(iIndividual) = var(data(isIndividual&hasCondition));
    dofIndividual(iIndividual) = length(data(isIndividual&hasCondition));     
    conditionIndividual(iIndividual) = unique(conditions(isIndividual&hasCondition));
  end

  mnCondition_unweighted(iCondition) = mean(mnIndividual);
  smCondition_unweighted(iCondition) = sem(mnIndividual);
end
% 
% % now follow the methods of moment estimator
% [conditionIndividualAll, mnIndividualAll, smIndividualAll, varIndividualAll, dofIndividualAll] = deal([]);
% for iCondition = 1:nConditions
%   
%   hasCondition = conditions==uniqueConditions(iCondition);  
%   uniqueIndividuals = unique(individuals(hasCondition));
%   nIndividuals      = length(uniqueIndividuals);
%   [mnIndividual,smIndividual, varIndividual,dofIndividual,conditionIndividual] = deal([]);
%   for iIndividual   = 1:nIndividuals
%     isIndividual = individuals==(uniqueIndividuals(iIndividual)); 
%     mnIndividual(iIndividual) = mean(data(isIndividual&hasCondition));
%     smIndividual(iIndividual) = sem(data(isIndividual&hasCondition));
%     varIndividual(iIndividual) = var(data(isIndividual&hasCondition));
%     dofIndividual(iIndividual) = length(data(isIndividual&hasCondition));     
%     conditionIndividual(iIndividual) = unique(conditions(isIndividual&hasCondition));
%   end
% 
%   mnIndividualAll  = [mnIndividualAll ; mnIndividual(:)];
%   smIndividualAll  = [smIndividualAll ; smIndividual(:)];
%   varIndividualAll = [varIndividualAll ; varIndividual(:)];
%   dofIndividualAll = [dofIndividualAll ; dofIndividual(:)];
%   conditionIndividualAll = [conditionIndividualAll ; conditionIndividual(:)];
% end
%   
% % randomize the conditions
% for iRandomization = 1:1000
%   perm = randperm(length(conditionIndividualAll));
%   conditionIndividualAll = conditionIndividualAll(perm);
%     
%   % now follow the methods of moment estimator
%   for iCondition = 1:nConditions
% 
%     hasCondition = conditionIndividualAll==iCondition;  
%     mnIndividual = mnIndividualAll(hasCondition);
%     smIndividual = smIndividualAll(hasCondition);
%     varIndividual = varIndividualAll(hasCondition);
%     dofIndividual = dofIndividualAll(hasCondition);
%     conditionIndividual = conditionIndividualAll(hasCondition);
%     
%      %%%%%%%%%%%%%%%%%%%%% compute the semi-weighted estimator
%     if strcmp(cfg.variance_estimator, 'pooled')
%       varPooled        = weightedmean(varIndividual,dofIndividual);
%       varIndividualEst = varPooled./dofIndividual;
%     else
%       varIndividualEst = varIndividual./dofIndividual;
%     end
% 
%     w = 1./varIndividualEst;
%     mu = weightedmean(mnIndividual,w);
% 
%     num = sum(w.*(mnIndividual-mu).^2) - (nIndividuals-1);
%     tausq(iCondition) = num ./ (sum(w) - sum(w.^2)./sum(w));
% 
%     Q(iCondition)    = sum(w.*(mnIndividual-mu).^2);
%     pval(iCondition) = compchi(Q(iCondition),nIndividuals-1); 
% 
%     if tausq(iCondition)<0, tausq(iCondition) = 0; end
% 
%     % compute the semiweighted estimator
%     weights      = 1./(varIndividualEst + tausq(iCondition));
%     mnCondition_semiweighted_rand(iCondition,iRandomization) = weightedmean(mnIndividual,weights);
%     smCondition_semiweighted_rand(iCondition,iRandomization) = sqrt(1./sum(weights));
%   end 
% end

% now follow the methods of moment estimator
for iCondition = 1:nConditions
  
  hasCondition = conditions==uniqueConditions(iCondition);  
  uniqueIndividuals = unique(individuals(hasCondition));
  nIndividuals      = length(uniqueIndividuals);
  [mnIndividual,smIndividual, varIndividual,dofIndividual,conditionIndividual] = deal([]);
  for iIndividual   = 1:nIndividuals
    isIndividual = individuals==(uniqueIndividuals(iIndividual)); 
    mnIndividual(iIndividual) = mean(data(isIndividual&hasCondition));
    smIndividual(iIndividual) = sem(data(isIndividual&hasCondition));
    varIndividual(iIndividual) = var(data(isIndividual&hasCondition));
    dofIndividual(iIndividual) = length(data(isIndividual&hasCondition));     
    conditionIndividual(iIndividual) = unique(conditions(isIndividual&hasCondition));
  end

  %%%%%%%%%%%%%%%%%%%%% compute the semi-weighted estimator
  if strcmp(cfg.variance_estimator, 'pooled')
    varPooled        = weightedmean(varIndividual,dofIndividual);
    varIndividualEst = varPooled./dofIndividual;
  else
    varIndividualEst = varIndividual./dofIndividual;
  end
  
    
  w = 1./varIndividualEst;
  worig = w;
  mu = weightedmean(mnIndividual,w);   
  Q(iCondition)    = sum(worig.*(mnIndividual-mu).^2);
  pval(iCondition) = compchi(Q(iCondition),nIndividuals-1); 

  if strcmp(cfg.semiweighted_adaptive,'yes')% max likelihood estimator
    for iIter = 1:cfg.semiweighted_adaptive_niter
      mu = weightedmean(mnIndividual,w);
      
      if iIter==1
        num = sum(w.*(mnIndividual-mu).^2) - (nIndividuals-1);        
        tausq(iCondition) = num ./ (sum(w) - sum(w.^2)./sum(w));
      else
        tausq(iCondition) = sum ( w.*( (mnIndividual-mu).^2 - varIndividualEst ) ) ./ sum(w);
      end
         
      if tausq(iCondition)<0, tausq(iCondition) = 0; end

      % compute the semiweighted estimator
      weights      = 1./(varIndividualEst + tausq(iCondition));
      mnCondition_semiweighted(iCondition) = weightedmean(mnIndividual,weights);
      smCondition_semiweighted(iCondition) = sqrt(1./sum(weights));
      nIndividuals_total(iCondition) = nIndividuals;
      nCells{iCondition} = dofIndividual;
      means{iCondition}  = mnIndividual;
      w = weights;
    end
  elseif strcmp(cfg.semiweighted_adaptive,'no')
    mu = weightedmean(mnIndividual,w);
    num = sum(w.*(mnIndividual-mu).^2) - (nIndividuals-1);

    if strcmp(cfg.tausq_est, 'weighted')
      tausq(iCondition) = num ./ (sum(w) - sum(w.^2)./sum(w));
    else
      tausq(iCondition) = var(mnIndividual) - (1./nIndividuals)*sum(w);
    end

    Q(iCondition)    = sum(w.*(mnIndividual-mu).^2);
    pval(iCondition) = compchi(Q(iCondition),nIndividuals-1); 

    if tausq(iCondition)<0, tausq(iCondition) = 0; end

    % compute the semiweighted estimator
    weights      = 1./(varIndividualEst + tausq(iCondition));
    mnCondition_semiweighted(iCondition) = weightedmean(mnIndividual,weights);
    smCondition_semiweighted(iCondition) = sqrt(1./sum(weights));
    nIndividuals_total(iCondition) = nIndividuals;
    nCells{iCondition} = dofIndividual;
    means{iCondition}  = mnIndividual;
  end
end 

% get some p values
% for iC = 1:nConditions-1
%    for jC = iC+1:nConditions
%      actD = abs(mnCondition_semiweighted(iC) - mnCondition_semiweighted(jC));
%      sD   = 
%    end
% end

% % determine the pvals of all the condition differences
% nConditions = length(mnCondition_semiweighted);
% pvals_all = NaN(nConditions,nConditions);
% for iC = 1:nConditions-1
%   for jC = iC+1:nConditions
%     actD = abs(mnCondition_semiweighted(iC) - mnCondition_semiweighted(jC));
%     randD = abs(mnCondition_semiweighted_rand(iC,:)-mnCondition_semiweighted_rand(jC,:));
%     randD_sort = sort(randD);
%     pvals_all(iC,jC) = 1-nearest(randD_sort, actD)/1000;
%   end
% end
    
%stat.pvals_all          = pvals_all;
stat.means_unweighted   = mnCondition_unweighted;
stat.sems_unweighted    = smCondition_unweighted;
stat.means_semiweighted = mnCondition_semiweighted;
stat.sems_semiweighted  = smCondition_semiweighted;
stat.means_weighted     = mnCondition_weighted;
stat.sems_weighted      = smCondition_weighted;
stat.means_semiweighted = mnCondition_semiweighted;
stat.sems_semiweighted  = smCondition_semiweighted;
stat.nIndividuals       = nIndividuals_total;
stat.nCells             = nCells;
stat.means              = means;
stat.tausq              = tausq;
stat.Q                  = Q;
stat.heterogeneity_pval = pval;


function mn = weightedmean(means,weights)

% for a vector
mn = sum(weights.*means) ./ sum(weights);

function sm = sem(x,dim)

if nargin < 2
	dim = find(size(x)~=1,1,'first');
	if isempty(dim)
		dim = 1; 
	end	  
end

sm = sqrt(var(x,[],dim)./size(x,dim));


function prob = compchi(Y,df)
    try
        prob = 1-chi2cdf(Y,df);%chance that null hypothesis (common direction) is right
%        disp('using the statistics toolbox')
    catch
        prob = chistat2(df,Y,inf);
    end

    