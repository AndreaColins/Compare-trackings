function MdlLinear=LDAlenght
%%%%this version of PCA analyse data in exploratory/touch periods, i.e,
%%%%times series may have different lenghts

%pole=char('smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
pole=char('air','smooth pole', 'carbon pole', 'black sandpaper', 'closed coil','open coil','Bamboo','Toothpick');
pole2=pole;
colors={'b','r','g','k','c','m','y'};
npole=size(pole,1);
nfiles=zeros(npole,1);
nresult=zeros(npole,1);
ms=50;%number of miliseconds
overlap=0;%miliseconds
poleaz=[];
poleel=[];
polekcor=[];
polekhor=[];
for i=1:npole
    %%%%%%%%%%%%%%%raw data
    filename=strcat('./tr4s/processed/FrameCorrected/result',pole(i,:),'.mat');
    %%%%%%%%%%%%%%%%%%%%this data is already frame corrected
    %filename=strcat('./tr4s/processed/bandpass/upto5/result',pole(i,:),'.mat');
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
        d=strcat('./tr4s/processed/touchdetection/',pole(i,:),'/');
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
%             result=centredk(endp,vtouch.touches,result);
%             if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
% 
%             result=result(startp:endp,:);
            
            %%%%%%%%%%%%for touch periods
             result=centredk(endp,vtouch.touches,result);
             result=touchselection(vtouch.touches,result);
             
            im=sum(imag(result(:)));
            if im>0
                display('Complex numbers found')
            end
            %%%%%%%%%%%for non touching periods
            %result=notouchselection(vtouch.touches,result);
            
            %%%%%%%%%for debugging
           % plot(result(:,4))
            %result=centredk(endp,vtouch.touches,result);
%              if (3488-endp)>=50
%                 endp=endp+50;
%             else
%                 endp=3488;
%             end
            %result=notouchselection(vtouch.touches,result);
            %result=touchselection(vtouch.touches,result);
                         %hold on
%                           result=result(startp:endp,:);
                         %plot(result(:,1),colors{i})
%                          hold off
%                         pause(5)
           
            
 
            
            nsamples=floor(size(result,1)/ms);
            nresult(i)=nresult(i)+nsamples;
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
    clear result
end
%%
azel=[tazimuth',televation'];
%%%%%%%%%Raw Data
size(azel)
plot(azel(1,:)','b')
hold on 
plot(azel(nresult(1)+1,:)','r')
%plot(azel(nresult(1)+nresult(2)+1,:)','g')
plot(azel(1:nresult(1),:)','b')
plot(azel(nresult(1)+1:nresult(1)+nresult(2),:)','r')
%plot(azel(nresult(1)+nresult(2)+1:end,:)','g')
hold off
%%

if strcmp(pole(1,1:3),'air')
    totalfiles=nfiles(1)*nresult(1)
    poleclass(1:totalfiles)={'air'};
    class1=azel(1:totalfiles,:);
else
    totalfiles=nresult(1)
    class1=azel(1:totalfiles,:);
    poleclass(1:totalfiles)={pole(1,:)};
end
for i=2:npole
    poleclass(totalfiles+1:totalfiles+nresult(i))={pole(i,:)};
    class2=azel(totalfiles+1:totalfiles+nresult(i),:);
totalfiles=totalfiles+nresult(i);

end
%%%%Apply PCA angles
MdlLinear = fitcdiscr(azel,poleclass);

MdlLinear.Coeffs(1,2).Linear



m1=mean(class1)';
m2=mean(class2)';
S_1=cov(class1);
S_2=cov(class2);
S_w=S_1+S_2;
S_w1=inv(S_w);

S_b=(m1-m2)*(m1-m2)';
[V,D]=eig(S_w1*S_b);

D=diag(D);
[Dsort,idx]=sort(D,'descend');


%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end
V=Vsort;
clear Vsort;

% figure
% scatter(class1(:,1),class1(:,2),'b')
% hold on 
% scatter(class2(:,1),class2(:,2),'r')

% f=@(x,y) V(1,2)*x+V(2,2)*y;
% h3= ezplot(f,[80 120 80 120]);
% h3.Color = 'g';
% hold off

class1p=V(:,1)'*class1';
class2p=V(:,1)'*class2';

figure
plot(class1p)
hold on
plot(class2p,'r')
hold off

figure
hi1=histogram(class1p)
hold on
hi2=histogram(class2p)
legend(pole(1,:),pole(2,:),'Location','best')
hold off
% size(SCORE)
% %%%%%%%%Plot results
% figure
% %%%graficar en los nuevos vectores
% subplot(2,2,1)
% scatter(SCORE(1:totalfiles,1),SCORE(1:totalfiles,2))
% for i=2:npole
% hold on
% scatter(SCORE(totalfiles+1:(totalfiles+nresult(i)),1),SCORE(totalfiles+1:(totalfiles+nresult(i)),2))   
% title('[Azimuth,elevation]')
% xlabel('First component')
% ylabel('Second component')
% totalfiles=totalfiles+nresult(i)
% end
% legend(pole2,'Location','best','Orientation','horizontal') 
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%plot Firsts components
% subplot(2,2,2)
% plot(COEFF(:,1),'r')
% hold on
% plot(COEFF(:,2),'b')
% legend('First Component','Second Component')
% xlabel('Time [miliseconds]')
% vexplained=EXPLAINED(1)+EXPLAINED(2);
% title(strcat('Variance explained=',num2str(vexplained)))
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%explained variance
% subplot(2,2,3)
% plot(cumsum(EXPLAINED),'*-')
% xlabel('Number of principal component')
% ylabel('Variance explained [%]')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%scatter of two principal components
% 
% subplot(2,2,4)
% scatter(COEFF(:,1),COEFF(:,2))
% xlabel('First Component')
% ylabel('Second Component')
% hold off


end
function [result]=centredk(endp,touches,result)
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
function [result]=touchselection(touches,result)
% plot(result(:,4))
% hold on
I=find(touches);
result=result(I,:);
% plot(result(:,4),'.r')
% pause(5)
% hold off
clear I
end