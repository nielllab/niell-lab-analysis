
badwv=wv(14,:); %enter value where waveform value discriminates bad events from good
bad=cells(badwv<.42,:)
b=find(badwv <.42)

%list of  2 5 8 10 16 18 
u=4;

%row org
cells(u,:) = [];
wn(u,:)=[];

psth(u,:)=[];
%layer(u,:)=[];
drift(u,:)=[];

%column org

% all_fit(:,u)=[];
% all_img(:,u)=[];
% all_test_img(:,u)=[];
% params(:,u)=[];

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