% Matlab codes for reading waveforms from TTank

clear;

ttx_ini;

figure
for ch = 1:16
    %N = invoke(TTX, 'ReadEventsV', 5000, Event_Name_Wave, 1, 0,  Select_Duration(1), Select_Duration(2), 'ALL')
    N = invoke(TTX, 'ReadEventsV', 5000, Event_Name_EED, 1, 0, 0,0, 'ALL')

    W = invoke(TTX, 'ParseEvV', 0, N);
    X=W(:);
    %plot(T_Wave(1: length(X)), X);

    subplot(4,4,ch);
    plot(X);
    axis([0 size(X/10,1) -100 100])
    axis off
end

    
    
TimeStamp = invoke(TTX, 'ParseEvInfoV', 0, N, 6);
%save (Save_File_Wave, 'TimeStamp', 'W');

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');
