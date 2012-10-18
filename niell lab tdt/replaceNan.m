function obs=replaceNan(obs,obs_name,field,value);
ind = find(strcmp(obs_name,field));
nanfind = isnan(obs(:,ind));
obs(nanfind,ind)=value;

