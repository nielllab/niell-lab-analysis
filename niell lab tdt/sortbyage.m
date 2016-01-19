function [sortedVals sortedError N]= sortbyage(vals,age,agelist,used);

for i = 1:length(agelist);
    i
if size(agelist,1)==1
    age_used = used & (age == agelist(i));
else
    age_used = used & (age >=agelist(1,i) & age<=agelist(2,i));
end
    N(i) = sum(age_used);
    sortedVals(i,:) = nanmean(vals(age_used,:),1);
    sortedError(i,:) = nanstd(vals(age_used,:),[],1)/sqrt(N(i));
    
end
