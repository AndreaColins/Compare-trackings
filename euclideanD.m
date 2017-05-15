function [rows,cols,val]=euclideanD(matrix)
distance=ones(size(matrix,2),size(matrix,2)).*1000;
for j=1:690
for i=j+1:size(matrix,2)
    distance(j,i)=norm(matrix(:,j)-matrix(:,i));
end
end
[rows,cols,val]=find(distance<=5);

end