function spots = addspot(sz,speeds,npersize,N,spots,nx,ny,cumprob,cumprobspeed);
%%% spots = [x y sz speed theta contrast]

side = ceil(rand(1)*4);
pos = ceil(rand(1)*nx);
r=rand(1);
spots(N,6) = round(rand)*2-1;
for i = 1:length(sz)
    
    if r<cumprob(i)
        spots(N,3)=(i);
        break
    end
end
r = rand(1);

for i = 1:length(speeds)
    
    if r<cumprobspeed(i)
        spots(N,4)=speeds(i);
        break
    end
end

if side ==1
    spots(N,1)=1;
    spots(N,2)=pos;
    spots(N,5) = pi*rand -pi/2;
elseif side ==2
    spots(N,1)=nx;
    spots(N,2)=pos;
    spots(N,5) = pi/2 +pi*rand;
elseif side ==3;
    spots(N,1)=pos;
    spots(N,2)=1;
    spots(N,5) = pi*rand;
elseif side==4;
    spots(N,1)=pos;
    spots(N,2)=ny;
    spots(N,5) = pi+pi*rand;
end

