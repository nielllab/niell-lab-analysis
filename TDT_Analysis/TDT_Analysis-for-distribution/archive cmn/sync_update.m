%%% updates for udp sync
%%% fill in to psychstim



%%% in setup  
stopudp = pnet('udpsocket',3787);
statusfile = fopen('statusfile.txt','w');
startTime = getsecs();

%%%during inter-stimulus break

if pnet(stopudp,'readpacket',20,'noblock')>0
    msg = pnet(stopudp,'read','char')
    if strcmp(msg,'stop')
        doneStim=1;
    end
end
elapsedTime = GetSecs-startTime;
fprintf(statusfile,'%d %f \n', iter,elapsedTime);

%%% at end
pnet('closeall')
fclose(statusfile)
    
    
    %%%  to send stop packet
udp = pnet('udpsocket',3787);
pnet(udp,'write','stop');
pnet(udp,'writepacket','mps-pc102',3787);
pnet('closeall')