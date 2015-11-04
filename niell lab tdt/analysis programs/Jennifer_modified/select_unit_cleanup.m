
badwv=wv(14,:); %enter value where waveform value discriminates bad events from good
bad=cells(badwv<.25,:)
find(badwv <.25)

%list of  11 5
u=5;

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