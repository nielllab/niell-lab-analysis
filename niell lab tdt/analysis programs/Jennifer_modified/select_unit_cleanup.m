
badwv=wv(13,:); %enter value where waveform value discriminates bad events from good
bad=cells(badwv<0.19,:)
find(badwv <0.19)

%list of 34 30 25 10 3 

u=7;

%row org
cells(u,:) = [];
wn(u,:)=[];

psth(u,:)=[];
layer(u,:)=[];
drift(u,:)=[];

%column org
trough_width(:,u) = [];
trough_depth(:,u) = [];
trough2peak(:,u)=[];
spikeT(:,u)=[];
peakchan(:,u)=[];
peak_height(:,u)=[];
nspikes(:,u)=[];
L_ratio(:,u)=[];
wv(:,u)=[];

save 'analysis_2.mat'