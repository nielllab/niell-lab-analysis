function plotDevData(vals, age, agelist, ageBins, used, label)

[resp resperr] = sortbyage(vals,age,agelist,used)
figure
errorbar(agelist,resp, resperr);
ylabel(label)

[resp resperr] = sortbyage(vals,age,ageBins,used)
figure
hold on
errorbar(mean(ageBins,1),resp, resperr);
bar(mean(ageBins,1),resp);
ylabel(label)