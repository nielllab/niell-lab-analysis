function show2dim_byclass(davec,clas,id)
% showing 2-dimensional data in X-Y space

nk=max(clas); 
divc=cell(1,nk);
pt=0;

for i=1:nk
    d=find(clas==i);
    divc{i}=d;
end

if nargin==3
  if id=='bk' 
      symb=['k.';'ks';'k+';'ko';'k*';'k^';'k.';'k+';'ko';'k*';'ks';'k.'];
  elseif id=='co' 
      symb=['b.';'rs';'mo';'g*';'k.';'c^';'ys';'b+';'mo';'g*';'k.';'c^'];
  elseif id=='nb'
      pt=1;
  else
      symb=[id;id;id;id;id;id;id;id;id;id;id;id;]; 
  end
else
  symb=['k.';'m.';'b.';'g.';'c.';'y.';'k+';'m+';'b+';'g+';'c+';'y+'];
end

%title('plot 2-dimensional data');
for i=1:nk
    hold on;
    d=divc{i};
    vh=max(davec);     
    vp=min(davec);
    if i==1 
      h1=vh(1); h2=vh(2); 
      p1=vp(1); p2=vp(2);       
    end
    if vh(1)>h1 
        h1=vh(1); 
    end
    if vh(2)<h2 
        h2=vh(2); 
    end
    if vp(1)>p1 
        p1=vp(1); 
    end
    if vp(2)<p2 
        p2=vp(2); 
    end
    
  if pt 
     plot(davec(d,1),davec(d,2),'y.');
     text(davec(d,1),davec(d,2),['\fontsize{10}' int2str(i)]);
   else
     plot(davec(d,1),davec(d,2),[symb(mod(i-1,12)+1,:)]); 
  end
end
box on;