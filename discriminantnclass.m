function discriminantnclass
%%%%inteto de hacer un LDA a mano

%%%two classes
class1=[3+(2*rand(50,1)-1),3+(2*rand(50,1)-1)];
class2=[-4+(2*rand(50,1)-1),6+(2*rand(50,1)-1)];
class3=[10+(2*rand(50,1)-1),10+(2*rand(50,1)-1)];


m1=mean(class1)';
m2=mean(class2)';
S_1=cov(class1);
S_2=cov(class2);
S_w=S_1+S_2;
S_w1=inv(S_w);
w=S_w1*(m1-m2);

S_b=(m1-m2)*(m1-m2)';
[V,D]=eig(S_w1*S_b)

D=diag(D);
[Dsort,idx]=sort(D,'descend');

%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end

V=Vsort;
figure
scatter(class1(:,1),class1(:,2),'b')
hold on 
scatter(class2(:,1),class2(:,2),'r')

f=@(x,y) V(1,2)*x+V(1,1)*y;
h3= ezplot(f,[-6 6 -6 6]);
h3.Color = 'g';
hold off

class1p=V(:,1)'*class1';
class2p=V(:,1)'*class2';

figure
plot(class1p)
hold on
plot(class2p,'r')
hold off

figure

hist(class1p)
hold on 
hist(class2p)
hold off

%%% k classes

scatter(class1(:,1),class1(:,2),'b')
hold on 
scatter(class2(:,1),class2(:,2),'r')
scatter(class3(:,1),class3(:,2),'g')

m1=mean(class1)';
m2=mean(class2)';
m3=mean(class3)';
mo=mean([class1;class2;class3])';
S_1=cov(class1);
S_2=cov(class2);
S_3=cov(class3);
S_w=S_1+S_2+S_3;
S_w1=inv(S_w);

S_b=((m1-mo)*(m1-mo)'+(m2-mo)*(m2-mo)'+(m3-mo)*(m3-mo)')/3;
[V,D]=eig(S_w1*S_b)


D=diag(D);
[Dsort,idx]=sort(D,'descend');

%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end

V=Vsort
class1p=V(:,1)'*class1';
class2p=V(:,1)'*class2';
class3p=V(:,1)'*class3';

t=-6:1:6;
wx=t.*V(1,1);
wy=t.*V(1,2);

plot(wx,wy)
figure
plot(class1p)
hold on
plot(class2p,'r')
plot(class3p,'g')
hold off

figure

histogram(class1p)
hold on 
histogram(class2p)
histogram(class3p)
hold off

LDA([class1;class2;class3],'scatter',{'class1','class2','class3'},[size(class1,1);size(class1,1);size(class1,1)])
end

function LDA(data,titlevar,classes,nresult)
totalfiles=nresult;
means=ones(size(data,2),size(nresult,1));
S_classes=means;
means(:,1)=mean(data(1:nresult(1)))';
S_classes=cov(data(1:nresult(1)));
for i=2:size(nresult,1)
    means(:,i)=mean(data(totalfiles+1:nresult(i)+totalfiles,:))';
    S_classes=S_classes+cov(data(totalfiles+1:nresult(i)+totalfiles,:));
    totalfiles=totalfiles+nresult(i);
end   

mo=mean(data)';
S_w=S_classes;
S_w1=inv(S_w);

%S_b=((m1-mo)*(m1-mo)'+(m2-mo)*(m2-mo)'+(m3-mo)*(m3-mo)')/3;
S_b=zeros(size(S_w));
for i=1:size(nresult,1)
    S_b=S_b+(means(:,i)-mo)*(means(:,i)-mo)';
    
end
S_b=S_b./size(nresult,1);
[V,D]=eig(S_w1*S_b);


D=diag(D);
[Dsort,idx]=sort(D,'descend');

%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end

V=Vsort;

classp=V(:,1)'*data';


figure
histogram(classp(1:nresult(1)))
totalfiles=nresult(1);
hold on
for i=2:size(nresult,1)
    histogram(classp(totalfiles+1:nresult(i)+totalfiles));
    totalfiles=totalfiles+nresult(i);
end   
title(titlevar)
legend(classes,'Location','best')
hold off
end