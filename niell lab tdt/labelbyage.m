function [vals bin]= sortbyage(vals,age,agelist,used);
bin = zeros(size(vals));
for i = 1:length(agelist);
    i
if size(agelist,1)==1
    age_used = used & (age == agelist(i));
else
    age_used = used & (age >=agelist(1,i) & age<=agelist(2,i));
end

bin(age_used) = i;
end

bin = bin(used);
vals = vals(used);
