figure
for i = 1:4
    lfp_raw(i,:) = lfp(i,:)./(1:size(lfp,2));
end
plot(f(f<58),lfp_raw(:,f<58)/500);
figure
semilogy(f(f<58 & f>2),lfp(:,f<58&f>2));
figure
semilogy(f(f<58 & f>2),lfp_raw(:,f<58&f>2)/max(max(lfp_raw(:,f<58&f>2))));
figure
loglog(f(f<58 & f>2),lfp_raw(:,f<58&f>2)/max(max(lfp_raw(:,f<58&f>2))));

figure
plot(f(f<58),lfp(:,f<58)/max(max((lfp(:,f<58)))));

