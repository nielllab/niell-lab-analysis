function [sortedVals sortedError N]= sortbyage(vals,age,agelist,used);
%%% binomial
for i = 1:length(agelist);
    i
if size(agelist,1)==1
    age_used = used & (age == agelist(i));
else
    age_used = used & (age >=agelist(1,i) & age<=agelist(2,i));
end

    N(i) = sum(age_used & ~isnan(vals(:,1)')); %%% added test for isnan 060217 cmn
    [mn ci] = binofit(sum(vals(age_used & ~isnan(vals(:,1)'))) ,sum(age_used & ~isnan(vals(:,1)'))');
    sortedVals(i,:) = mn;
    sortedError(i,1) = ci(1);
    sortedError(i,2)= ci(2);
    
end
