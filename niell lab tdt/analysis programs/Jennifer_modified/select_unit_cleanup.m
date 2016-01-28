
badwv=wv(14,:); %enter value where waveform value discriminates bad events from good
bad=cells(badwv<.29,:)
b=find(badwv <.29)

%list of  5 4 3 2 1 
u=1;

%row org
cells(u,:) = [];
%wn(u,:)=[];

%psth(u,:)=[];
%layer(u,:)=[];
%drift(u,:)=[];

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