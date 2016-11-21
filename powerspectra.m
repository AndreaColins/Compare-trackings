function powerspectra(vectors)
n=size(vectors,2);
t=size(vectors,1);
fftvector=zeros(t,n);
for i=1:n
    fftvector(:,i)=(abs(fft(vectors(:,i)))).^2/t;
end
semilogy(fftvector)
hold on 
legend('Azimuth','Elevation','Kappa Coronal','Kappa Horizontal')
hold off
clear n t fftvector 
end