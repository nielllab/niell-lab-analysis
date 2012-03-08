for f = 1:300;
    imshow(moviedata(:,:,f));
    mov(f) = getframe(gcf);
end
movie2avi(mov,'noisemovie.avi','Fps',30);
