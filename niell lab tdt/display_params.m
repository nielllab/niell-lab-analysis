function display_params

[afname, apname] = uigetfile('*.mat','analysis data');
prompt = {'monitor offset deg (midline=0) : ','monitor height (cm) : '};
num_lines = 1;
def = {'45','10'};
answer = inputdlg(prompt,'display parameters',num_lines,def);
displayOffset = str2num(answer{1})
displayHeight = str2num(answer{2})
save(fullfile(apname,afname),'displayOffset','displayHeight','-append');
