function powerspectra(vectors)
n=size(vectors,2)
t=size(vectors,1)
fftvector=zeros(t,n);
for i=1:n
    fftvector(:,i)=(abs(fft(vectors(:,i)))).^2/t;
end
loglog(fftvector)
hold on 
legend('1','2','3','4')
hold off
clear n t fftvector 
end