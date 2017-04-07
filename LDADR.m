function [classp,azel]=LDADR
%%%%this version of LDA analyse data in exploratory/touch periods, i.e, and
%%%%is used to dimensionality reduction in times series may have different lenghts

pole=char('air','Smooth pole', 'Carbon Pole', 'Black Sandpaper', 'Closed coil','Open coil','Cardboard','Bamboo','Toothpick','Wood');
%pole=char('air','Smooth pole');
%pole=char('air','smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms=20;%number of miliseconds
overlap=0;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
resultp=[];
for i=1:npole
    %%%%%%%%%%%%%%%raw data random selection
    filename=strcat('./videosselectedtextures/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%raw data
    %filename=strcat('./tr4s/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./tr4s/processed/bandpass/from20/result',pole(i,:),'.mat');
    if strcmp(pole(i,1:3),'air')
        v=load(filename,'-mat');
        result=v.result;
        nfiles(i)=size(result,3);
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
        
    else 
        d=strcat('./videosselectedtextures/processed/touchdetection/',pole(i,:),'/');
        ff = dir([d '*.mat']);
        nfiles(i)=size(ff,1);
        string_list=cell(nfiles(i),1);
        for j=1:nfiles(i)
            string_list{j}=strcat(d,ff(j).name);
            touchfile=strcat(d,ff(j).name)
            v=load(filename,'-mat');
            result=v.result;
            vtouch=load(touchfile,'-mat');
            idx=[1:3488];
            idx2=circshift(idx,size(idx,1)-vtouch.start_frame-1,2);
            result=result(idx,:,j);
            vtouch.touches=vtouch.touches(idx2);
            %%%%%%%select touch periods
            startp=find(vtouch.touches,1,'first');
            endp=find(vtouch.touches,1,'last');

            %%%%%%%%%%%%%%%for exploratory periods
%             [result,touchperiods]=centredk(endp,vtouch.touches,result);
%             if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
% 
%             result=result(startp:endp,:);
%             
            %%%%%%%%%%%%for touch periods
             [result,touchperiods]=centredk(endp,vtouch.touches,result);
             result=touchselection(vtouch.touches,result,touchperiods);
             
            im=sum(imag(result(:)));
            if im>0
                display('Complex numbers found')
            end
            %%%%%%%%%%%for non touching periods
            %result=notouchselection(vtouch.touches,result);
            
            %%%%%%%%%for debugging
%            plot(result(:,4))
%             result=centredk(endp,vtouch.touches,result);
%              if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
            %result=notouchselection(vtouch.touches,result);
            %result=touchselection(vtouch.touches,result);
%                          hold on
%                          result=result(startp:endp,:);
%                          plot(result(:,4),colors{i})
%                          hold off
%                          pause(5)
          
           
            %resultp=[resultp; result];
%             plot(resultp(:,1))
%             pause(2)
%             clear result
        %end
        %result=resultp;
            nsamples=floor(size(result,1)/ms);
            nresult(i)=nsamples;
            result((nsamples*ms+1):end,:,:)=[];
    
            poleaz2=reshape(squeeze(result(:,1,:)),nsamples.*ms,1);
            poleaz2=reshape(poleaz2,ms,nsamples);
          
            
            poleel2=reshape(squeeze(result(:,2,:)),nsamples.*ms,1);
            poleel2=reshape(poleel2,ms,nsamples);
            
            polekcor2=reshape(squeeze(result(:,3,:)),nsamples.*ms,1);
            polekcor2=reshape(polekcor2,ms,nsamples);
            
            polekhor2=reshape(squeeze(result(:,4,:)),nsamples.*ms,1);
            polekhor2=reshape(polekhor2,ms,nsamples);
            
            poleaz=[poleaz,poleaz2];
            poleel=[poleel,poleel2];
            polekcor=[polekcor,polekcor2];
            polekhor=[polekhor,polekhor2];
     
        end 
        
    end
   
   
    
    %%%%%%%%%%%%%%%%%%Add new condition
    if i==1
        tazimuth=poleaz;
        televation=poleel;
        tkcoronal=polekcor;
        tkhorizontal=polekhor;
        
    else
        tazimuth=[tazimuth,poleaz];
        televation=[televation,poleel];
        tkcoronal=[tkcoronal,polekcor];
        tkhorizontal=[tkhorizontal,polekhor];
    end
    
    poleaz=[];
    poleel=[];
    polekcor=[];
    polekhor=[];
    resultp=[];
    clear result 
end
%%
%%%[Azimuth,Elevation]
azel=[zscore(tazimuth'),zscore(televation')];

totalfiles=nfiles(1)*nresult(1)
nresult(1)=totalfiles;

poleclass(1:totalfiles)={'air'};
class1=azel(1:totalfiles,:);
poleclass(totalfiles+1:size(azel,1))={'pole'};
class2=azel(totalfiles+1:end,:);
%%plot
%%%%%%%%%Raw Data
size(azel)

plot(azel(1,:)','b')
hold on
plot(azel(nresult(1)+1,:)','r')
plot(azel(1:nresult(1),:)','b')
plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
hold off
%%%%Apply LDA angles
%MdlLinear = fitcdiscr(azel,poleclass);

classp=LDA(class1,class2,'[Azimuth,Elevation]');

% %%%[Azimuth]
% azel=[tazimuth'];
% totalfiles=nfiles(1)*nresult(1);
% class1=azel(1:totalfiles,:);
% class2=azel(totalfiles+1:end,:);
% %%plot
% %%%%%%%%%Raw Data
% figure
% plot(azel(1,:)','b')
% hold on
% plot(azel(nresult(1)+1,:)','r')
% plot(azel(1:nresult(1),:)','b')
% plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
% hold off
% %%%%Apply LDA angles
% %MdlLinear = fitcdiscr(azel,poleclass);
% 
% LDA(class1,class2,'[Azimuth]')
% 
% %%%[Azimuth]
% azel=[televation'];
% totalfiles=nfiles(1)*nresult(1)
% class1=azel(1:totalfiles,:);
% class2=azel(totalfiles+1:end,:);
% %%plot
% %%%%%%%%%Raw Data
% figure
% plot(azel(1,:)','b')
% hold on
% plot(azel(nresult(1)+1,:)','r')
% plot(azel(1:nresult(1),:)','b')
% plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
% hold off
% %%%%Apply LDA angles
% %MdlLinear = fitcdiscr(azel,poleclass);
% 
% LDA(class1,class2,'[Elevation]')

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

function [classp]=LDA(class1,class2,titlevar)
m1=mean(class1)';
m2=mean(class2)';
S_1=cov(class1);
S_2=cov(class2);
S_w=S_1+S_2;
%S_w1=inv(S_w);

S_b=(m1-m2)*(m1-m2)';
[V,D]=eig(S_w\S_b);

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
title(titlevar)
legend('air','pole','Location','best')
hold off
subplot(2,1,2)
plot(V(:,1),'-*')
hold on 
xlabel('Time [ms]')
ylabel('First LDA component value')
hold off
classp=[class1p';class2p'];
end