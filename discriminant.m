function discriminant
%%%%inteto de hacer un LDA a mano

%%%two classes
class1=[3+(2*rand(50,1)-1),3+(2*rand(50,1)-1)];
class2=[-4+(2*rand(50,1)-1),4+(2*rand(50,1)-1)];

% scatter(class1(:,1),class1(:,2),'b')
% hold on 
% scatter(class2(:,1),class2(:,2),'r')


m1=mean(class1)';
m2=mean(class2)';
S_1=cov(class1);
S_2=cov(class2);
S_w=S_1+S_2;
S_w1=inv(S_w);
%w=S_w1*(m1-m2);

% f=@(x,y) w(2)*x+w(1)*y;
% f2=@(x,y) w(1)*x+w(2)*y;
% h2 = ezplot(f,[-6 6 -6 6]);
% h2.Color = 'r';
% h2.LineWidth = 2;
% h3= ezplot(f2,[-6 6 -6 6]);
% h3.Color = 'g';
% hold off

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
clear Vsort;

figure
scatter(class1(:,1),class1(:,2),'b')
hold on 
scatter(class2(:,1),class2(:,2),'r')

f=@(x,y) V(1,2)*x+V(2,2)*y;
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
%%% k classes

end