function LDAanalysis
% [classp,azel]=LDADRtextures;
% % [minp,imin]=min(classp);
% % [maxp,imax]=max(classp);
% 
% [sortedp,I]=sort(classp,'ascend');
% %meansample=mean(azel(:,size(azel,2)/2+1:end),2)-mean(azel(:,1:size(azel,2)/2),2);
% [env,~]=envelope(azel);
% meansample=mean(azel,2);
% [sortedm,Im]=sort(meansample,'ascend');
% figure
% subplot(3,1,1)
% azel=azel';
% plot(azel(:,I(end-5:end)))
% title('maximum score of texture')
% %axis([0 size(azel,1) -3 4])
% subplot(3,1,2)
% middle=floor(size(I,2)/2)-2
% plot(azel(:,I(middle-2:middle+2)))
% %axis([0 size(azel,1) -3 4])
% subplot(3,1,3)
% 
% plot(azel(:,I(1:5)))
% title('minimum score of whisking')
% %axis([0 size(azel,1) -3 4])
% 
% figure 
% scatter(I,Im)
tic
final=10;
ms=10;

obs=[ms:10:final]
nair=zeros(size(obs))';
npole=nair;
nrep=5;
misrate=zeros(nrep,size(obs,2));


j=1;
for i=obs
    fileID = fopen('classifier.txt','w');
    setGlobalID(fileID);
    i
    for rep=1:nrep
        
        [classp,azel,confusion,nair(j,1),npole(j,1)]=LDADRinprogress(i);
        misrate(rep,j)=(confusion(1,2)+confusion(2,1))/sum(sum(confusion));
        
    end
    fclose(fileID);
    fileID2 = fopen('classifier.txt','r');
    classifier = fscanf(fileID2,'%f');
    size(classifier)
    classifier=reshape(classifier,i,nrep*10)';
    errorbar(mean(classifier(:,1:end)),std(classifier(:,1:end)))
    hold on
    xlabel('Time [ms]')
    ylabel('Eigenvector')
    pause(10)
    fclose(fileID2);
    j=j+1;
    clear classifier
end

[minmis,xmin]=min(misrate);
minms=xmin+ms-1;
subplot(1,2,1)
misrate
errorbar(2*obs,mean(misrate,1),std(misrate,1))
hold on
xlabel('Windows size [ms]')
ylabel('Misclassification rate')
hold off
subplot(1,2,2)
ybar=[nair,npole]
bar(2*obs,ybar)
xlabel('Windows size [ms]');
ylabel('Number of samples')
legend('air','pole')
toc
end
function setGlobalID(val)
global x
x = val;
end
function r = getGlobalID
global x
r = x;
end