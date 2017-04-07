function reconstrucction(result)
N=3488;

n2=100;
t= [0:N-1]'/(N); 
T=t*3.488;
ft=fft(result);
a=-1/(N/2)*imag(ft);
a=a(1:n2);
b=1/(N/2)*real(ft);
b=b(1:n2);
rec=zeros(size(result))+b(1)/2;
for i=2:size(a,1)
    rec=rec+a(i)*sin(2*pi*(i-1)*t)+b(i)*cos(2*pi*(i-1)*t);
end
plot(T,rec,'*')
hold on
plot(T,result,'r')
%axis([0 1 -1.5 1.5])
legend('Original signal','Reconstructed signal')
hold off
end