function hw=fwhh(pole)

filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'.mat')
v=load(filename,'-mat');
nfiles=size(v.autoazimuth,2);
N0=round(size(v.autoazimuth(:,1),1)/2);
hw=zeros(nfiles,4);
auto=ones(nfiles,4);
figure

for i=1:nfiles
auto(i,1)=find(v.autoazimuth(:,i)>=0.5,1);
auto(i,2)=find(v.autoelevation(:,i)>=0.5,1);
auto(i,3)=find(v.autokcoronal(:,i)>=0.5,1);
auto(i,4)=find(v.autokhorizontal(:,i)>=0.5,1);
hw(i,:)=(N0*ones(1,4)-auto(i,:))*2;

plot(v.autoazimuth(:,i))
hold on
x2=auto(i,:):1:(auto(i,:)+hw(i,1));
y2=ones(size(x2))*0.5;
plot(x2,y2)
axis([3000 4000 0 1])
hold off
pause(0.5)


end
return