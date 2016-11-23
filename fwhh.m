function hw=fwhh(pole)

filename=strcat('./tr4s/tr4s/',pole,'/autocov',pole,'.mat')
v=load(filename,'-mat');
nfiles=size(v.autoazimuth,2);
N0=round(size(v.autoazimuth(:,1),1)/2);
hw=zeros(nfiles,4);
% figure
% plot(v.autoazimuth(:,1:4))
auto=ones(nfiles,4);

for i=1:nfiles
auto(i,1)=find(v.autoazimuth(:,i)>=0.5,1);
auto(i,2)=find(v.autoelevation(:,i)>=0.5,1);
auto(i,3)=find(v.autokcoronal(:,i)>=0.5,1);
auto(i,4)=find(v.autokhorizontal(:,i)>=0.5,1);
i


hw(i,:)=(N0*ones(1,4)-auto(i,:))*2;
end
return