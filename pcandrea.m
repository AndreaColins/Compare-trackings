function [Vsort,y,Dsort,varcorr,explained]=pcandrea(vectors)
vectors=zscore(vectors);
%cov calculate the covariance of every column 
c=cov(vectors);
%calculate eigenvectors of c (every column is an eigenvector)
% and eigenvalues D (matrix where every element in the diagonal is an
% eigenvalue)
[V,D]=eig(c);
%sort the eigenvalues in descendent order
D=diag(D);
[Dsort,idx]=sort(D,'descend');


%sort the corresponding eigenvectors
Vsort=zeros(size(V));

for i=1:size(V,1)
    Vsort(:,i)=V(:,idx(i));
end
explained=Dsort/sum(Dsort)*100;
%plot the proportion of variance explained by the eigen values to decide how many principal components choose
% plot(Dsort/sum(Dsort)*100)
% hold on
% xlabel('Number of component')
% ylabel('Percentage of variance explained')
% hold off

%calculate principal component scores
y=vectors*Vsort;

%correlate principal components with every variable (rows are variables and
%columns are number of component)
varcorr=corr(vectors,y);
end