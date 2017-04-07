%%%%hilbert transform

x=0:0.01:1;
y=sin(4*3.14*x)+0.5*sin(8*3.14*x);

figure
plot(x,y)
hold on 
h1=hilbert(y);
amplitude=abs(h1);

 plot(x,amplitude','r')
 
 figure
 histogram(amplitude,60)
 title('')
% 
% angle=atan(real(h1)./imag(h1));
% figure
% plot(angle)
