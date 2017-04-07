function [classp,azel,confusion,nair,npole]=LDADRinprogress
%%%%this version of LDA analyse data in exploratory/touch periods, i.e, and
%%%%is used to dimensionality reduction in times series may have different lenghts

pole=char('air','smooth pole','carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
%pole=char('smooth pole','carbon pole');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms1=10;%number of miliseconds
overlap=0;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
tazimuth=[];
televation=[];
tkcoronal=[];
tkhorizontal=[];
for i=1:npole
    ms=ms1;%number of miliseconds
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./videosselectedtextures/processed/bandpass/from20/result',pole(i,:),'.mat');
    if strcmp(pole(i,1:3),'air')
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
        ms=ms-overlap;
        nresult(i)=floor(size(result,1)/ms);
        result((nresult(i)*ms+1):end,:,:)=[];
        
        
        poleaz=reshape(squeeze(result(:,1,:)),nfiles(i).*nresult(i).*ms,1);
        poleaz=reshape(poleaz,ms,nresult(i).*nfiles(i));
        
        poleel=reshape(squeeze(result(:,2,:)),nfiles(i).*nresult(i).*ms,1);
        poleel=reshape(poleel,ms,nresult(i).*nfiles(i));
        
        polekcor=reshape(squeeze(result(:,3,:)),nfiles(i).*nresult(i).*ms,1);
        polekcor=reshape(polekcor,ms,nresult(i).*nfiles(i));
        
        polekhor=reshape(squeeze(result(:,4,:)),nfiles(i).*nresult(i).*ms,1);
        polekhor=reshape(polekhor,ms,nresult(i).*nfiles(i));
        
        if overlap~=0
            for j=1:nresult(i)*nfiles(i)-1
                poleaz2(:,j)=[poleaz(:,j);poleaz(1:overlap,j+1)];
                poleel2(:,j)=[poleel(:,j);poleel(1:overlap,j+1)];
                polekcor2(:,j)=[polekcor(:,j);polekcor(1:overlap,j+1)];
                polekhor2(:,j)=[polekhor(:,j);polekhor(1:overlap,j+1)];
            end
           
            poleaz= poleaz2;
            poleel=poleel2;
            polekcor=polekcor2;
            polekhor=polekhor2;
        end
    else 
        d=strcat('./videosselectedtextures/processed/touchdetection/',pole(i,:),'/');
        ff = dir([d '*.mat']);
        nfiles(i)=size(ff,1);
        string_list=cell(nfiles(i),1);
        for j=1:nfiles(i)
            string_list{j}=strcat(d,ff(j).name);
            touchfile=strcat(d,ff(j).name);
            v=load(filename,'-mat');
            result=v.result;
            vtouch=load(touchfile,'-mat');
            idx=[1:3488];
            idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
            %idx2=circshift(idx,-(size(idx,1)-vtouch.start_frame-1));
            result=result(idx,:,j);
            vtouch.touches=vtouch.touches(idx2);
            %%%%%%%select touch periods
            startp=find(vtouch.touches,1,'first');
            endp=find(vtouch.touches,1,'last');

            %%%%%%%%%%%%%%%for exploratory periods
%             [result2,touchp]=centredk(endp,vtouch.touches,result);
%             if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
%             touchp=[startp endp];
            %result=result2(startp:endp,:);
%             
            %%%%%%%%%%%%for touch periods
             [result2,touchp]=centredk(endp,vtouch.touches,result);
             %discard small touch periods
             counter=1;
             touchperiods=[];
             for t=1:size(touchp,1)
                 if (touchp(t,2)-touchp(t,1))>=ms
                     touchperiods(counter,:)=touchp(t,:);
                     counter=counter+1;
                 end
             end
           
             clear touchp 
             if ~isempty(touchperiods)
            for t=1:size(touchperiods,1)
                result=result2(touchperiods(t,1):touchperiods(t,2),:);

            nsamples=floor(size(result,1)/(ms-overlap));
            if mod(size(result,1),(ms-overlap))<overlap
                nsamples=nsamples-1;
            end
            nresult(i)=nsamples+nresult(i);
            result((nsamples*(ms-overlap)+overlap+1):end,:)=[];
            
             if overlap==0

                poleaz2=reshape(result(:,1),ms,nsamples);
                poleel2=reshape(result(:,2),ms,nsamples);
                polekcor2=reshape(result(:,3),ms,nsamples);
                polekhor2=reshape(result(:,4),ms,nsamples);
             else
                
                for j=0:nsamples-1
                    poleaz2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,1);
%                     plot(j*(ms-overlap)+1:1:(j+1)*ms-j*overlap,poleaz2(:,j+1))
%                     pause(2)
                    poleel2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,2);
                    polekcor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,3);
                    polekhor2(:,j+1)=result(j*(ms-overlap)+1:(j+1)*ms-j*overlap,4);
                end
                
            end
    
            poleaz=[poleaz, poleaz2];

            poleel=[poleel,poleel2];
            
            polekcor=[polekcor,polekcor2];
            
            polekhor=[polekhor, polekhor2];
            
            clear poleaz2 poleel2 polekcor2 polekhor2 result 
            end
             end
            clear touchperiods
        end 
        
    end
   
   
    
    %%%%%%%%%%%%%%%%%%Add new condition

        tazimuth=[tazimuth,poleaz];
        televation=[televation,poleel];
        tkcoronal=[tkcoronal,polekcor];
        tkhorizontal=[tkhorizontal,polekhor];

    
    poleaz=[];
    poleel=[];
    polekcor=[];
    polekhor=[];
    clear result poleaz2 poleel2 polekcor2 polekhor2 result2
end
%%
%%%[Azimuth,Elevation]
azel=[tazimuth',televation'];

%size(azel)
totalfiles=nresult(1)*nfiles(1);
nresult(1)=totalfiles;
poleclass(1:totalfiles)={'air'};
class1=azel(1:totalfiles,:);
poleclass(totalfiles+1:size(azel,1))={'pole'};
class2=azel(totalfiles+1:end,:);
%%plot
%%%%%%%%%Raw Data
% size(azel)
% 
% plot(azel(1,:)','b')
% hold on
% plot(azel(nresult(1)+1,:)','r')
% %plot(azel(1:nresult(1),:)','b')
% %plot(azel(nresult(1)+1:end,:)','r')
% hold off
%%%%Apply LDA angles
%MdlLinear = fitcdiscr(azel,poleclass);
nair=size(class1,1);
npole=size(class2,1);
poleclass=poleclass';
classp=LDA(azel,poleclass);
%confusion=crossval2(class1,class2,'[Azimuth,Elevation]',4);
%misrate=(confusion(1,2)+confusion(2,1))/sum(sum(confusion));

%%%%%%%%%%%debugging
% test=azel([totalfiles-2:totalfiles+2,5800:5805],:);
% azel([totalfiles-2:totalfiles+2,5800:5805],:)=[];
% poleclasstest=poleclass([totalfiles-2:totalfiles+2,5800:5805],:)
% poleclass([totalfiles-2:totalfiles+2,5800:5805],:)=[];
%confusion=crossval3(azel,poleclass,test,poleclasstest)


% cp = cvpartition(poleclass,'k',10);
% f=@(xtr,ytr,xtest,ytest)crossval3(xtr,ytr,xtest,ytest);
% confusion = crossval(f,azel,poleclass,'partition',cp);
% confusion=reshape(sum(confusion),2,2);
% misrate=(confusion(1,2)+confusion(2,1))/sum(sum(confusion));
% classp=0;

 


end
function [result,touchperiods]=centredk(endp,touches,result)
i=1;
touchescopy=touches;
touchcount=0;
while i<endp
    ini=find(touchescopy,1,'first');
    touchescopy(1:ini)=1;
    fin=find(~touchescopy,1,'first')-1;
     if isempty(fin)
        fin=3488;
    end
    touchescopy(1:fin)=0;
    hold on
    touchperiods(touchcount+1,1:2)=[ini,fin];
    touchcount=touchcount+1;
    i=fin;
end
%touchperiods

kappacoronal=result(:,3,1);
kappahorizontal=result(:,4,1);
kappacoronal(~touches)=0;
kappahorizontal(~touches)=0;
for i=1:touchcount
    if i==1
        if touchperiods(1,1)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(1:touchperiods(i,1)-1,4));
        end
    else
        if touchperiods(i,1)-touchperiods(i-1,2)>10
            k0coronal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i,1)-10:touchperiods(i,1)-1,4));
        else
            k0coronal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,3));
            k0horizontal=median(result(touchperiods(i-1,2)+1:touchperiods(i,1)-1,4));
        end
    end
        kappacoronal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),3)-k0coronal;
    kappahorizontal(touchperiods(i,1):touchperiods(i,2))=result(touchperiods(i,1):touchperiods(i,2),4)-k0horizontal;
   
end
result(:,3)=kappacoronal;
result(:,4)=kappahorizontal;
clear touchescopy ini fin
end
function [result]=notouchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches,1,'first');
I2=find(touches,1,'last');
if I2==3488
result=result(1:I,:);
else
result=[result(1:I,:);result(I2:end,:)];
end
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I I2
end
function [result]=touchselection(touches,result,touchperiods)
% plot(result(:,4))
% hold on
I=find(touches);
result=result(I,:);
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I
end

function [classp,classifier]=LDA(class,poleclass)
class1=[];
class2=[];
for i=1:size(class,1)
    if strcmp(poleclass(i,:),'air')
        class1=[class1;class(i,:)];
    else
        class2=[class2;class(i,:)];
    end
end
m1=mean(class1)';
m2=mean(class2)';
S_1=cov(class1);
S_2=cov(class2);
S_w=S_1+S_2;
%S_w1=inv(S_w);
S_b=(m1-m2)*(m1-m2)';
A=S_w\S_b;

[V,D]=eig(A);
sample=A*V-V*D;
sample(1,1);
D=diag(D);
[Dsort,idx]=sort(D,'descend');


%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end
V=Vsort;
clear Vsort;
V(:,1)=real(V(:,1));

class1p=V(:,1)'*class1';
class2p=V(:,1)'*class2';

% figure
% plot(class1p)
% hold on
% plot(class2p,'r')
% hold off

figure
subplot(2,1,1)
[~,edges]=histcounts([class1p';class2p'],20);
hold on
histogram(class1p,edges,'Normalization','pdf');
histogram(class2p,edges,'Normalization','pdf');
% [~,edges]=hist([class1p';class2p'],20);
% h1=hist(class1p,edges);
% h2=hist(class2p,edges);
% bar(h2,'r')
% hold on
% bar(h1,'b')

%title(titlevar)
legend('pole','air','Location','best')
hold off
subplot(2,1,2)
plot(V(:,1)./sqrt(sum(V(:,1).^2)),'-*')
%plot([class1p';class2p'],'-*')
hold on 
xlabel('Time [ms]')
ylabel('First LDA component value')
hold off
% subplot(3,1,3)
% [~,edges]=histcounts(mean([class1;class2],2),20);
% hold on
% histogram(mean(class1,2),edges,'Normalization','pdf');
% histogram(mean(class2,2),edges,'Normalization','pdf');
% [~,edges]=hist([class1p';class2p'],20);
% h1=hist(class1p,edges);
% h2=hist(class2p,edges);
% bar(h2,'r')
% hold on
% bar(h1,'b')

%title(titlevar)
%legend('pole','air','Location','best')
%hold off
classp=[class1p';class2p'];
classifier=real(V(:,1));
end

function confusion=crossval2(class1,class2,titlevar,kfold)
confusion=zeros(2,2);
nc1=floor(size(class1,1)/kfold);
nc2=floor(size(class2,1)/kfold);
clf=zeros(size(class1,2),1);
for i=1:kfold
    %%%%%partition
    class1train=class1;
    class2train=class2;
    class1test=class1(nc1*(i-1)+1:nc1*i,:);
    class1train(nc1*(i-1)+1:nc1*i,:)=[];
    class2test=class2(nc2*(i-1)+1:nc2*i,:);
    class2train(nc2*(i-1)+1:nc2*i,:)=[];
    %%%%%%training
    [classp,classifier]=LDA(class1train,class2train,titlevar);

    discr=(mean(classp(1:size(class1train)))+mean(classp(size(class1train)+1:end)))/2;
    
    %test
    class1ptest=classifier'*class1test';
    class2ptest=classifier'*class2test';
    
    if mean(classp(1:size(class1train)))>mean(classp(size(class1train)+1:end))%%%%%if class1 > class2
       %disp('air>pole')
        confusion(1,1)=size(find(class1ptest-discr>0)',1)+confusion(1,1);%%%%classified as class1 and it is class1
        confusion(1,2)=size(find(class1ptest-discr<0)',1)+confusion(1,2);%%%%classified as class2 and it is class1
        confusion(2,2)=size(find(class2ptest-discr<0)',1)+confusion(2,2);%%%%classified as class2 and it is class2
        confusion(2,1)=size(find(class2ptest-discr>0)',1)+confusion(2,1);%%%%classified as class2 and it is class1
    else
       %disp('pole>air')
        confusion(1,1)=size(find(class1ptest-discr<0)',1)+confusion(1,1);%%%%classified as class1 and it is class1
        confusion(1,2)=size(find(class1ptest-discr>0)',1)+confusion(1,2);%%%%classified as class2 and it is class1
        confusion(2,2)=size(find(class2ptest-discr>0)',1)+confusion(2,2);%%%%classified as class2 and it is class2
        confusion(2,1)=size(find(class2ptest-discr<0)',1)+confusion(2,1);%%%%classified as class2 and it is class2
    end
clf=clf+classifier;
    
end
% figure
% plot(clf/kfold)
end
function confusion=crossval3(classtrain,poleclasstr,classtest,poleclasstest)
confusion=zeros(2,2);

    %%%%%%training
    [classp,classifier]=LDA(classtrain,poleclasstr);
    fileID=getGlobalID;
   
    fprintf(fileID,'%f\n',classifier);
    I=find(strcmp(poleclasstr,'air'));
    meanclass1=mean(classp(I));
    I=find(~strcmp(poleclasstr,'air'));
    meanclass2=mean(classp(I));
    discr=mean([meanclass1,meanclass2]);
   
    %test
    classptest=classifier'*classtest';
    poleclassptest=poleclasstest;
    poleclassptest(1:end,:)={'pole'};
    I=find(classptest-discr>0);
    
    I2=find(classptest-discr<0);
    
    if meanclass1>meanclass2%%%%%if class1 > class2
        poleclassptest(I)={'air'};
        confusion(1,1)=size(find(strcmp(poleclassptest(I),poleclasstest(I))==1),1);%%%%classified as class1 and it is class1
        confusion(2,1)=size(find(strcmp(poleclassptest(I),poleclasstest(I))==0),1);%%%%classified as class2 and it is class1
        confusion(1,2)=size(find(strcmp(poleclassptest(I2),poleclasstest(I2))==0),1);%%%%classified as class2 and it is class2
        confusion(2,2)=size(find(strcmp(poleclassptest(I2),poleclasstest(I2))==1),1);%%%%classified as class2 and it is class1
        
    else
        poleclassptest(I2)={'air'};
        %strcmp(poleclassptest(I),poleclasstest(I))
         confusion(1,1)=size(find(strcmp(poleclassptest(I2),poleclasstest(I2))==1),1);%%%%classified as class1 and it is class1
        confusion(2,1)=size(find(strcmp(poleclassptest(I2),poleclasstest(I2))==0),1);%%%%classified as class2 and it is class1
        confusion(1,2)=size(find(strcmp(poleclassptest(I),poleclasstest(I))==0),1);%%%%classified as class2 and it is class2
        confusion(2,2)=size(find(strcmp(poleclassptest(I),poleclasstest(I))==1),1);%%%%classified as class2 and it is class1
       
    end


end
function r = getGlobalID
global x
r = x;
end

