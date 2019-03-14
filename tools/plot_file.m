clear
fID = fopen('pendulum.txt','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fID, formatSpec, sizeA);
A = A'
plot(A(:,1),A(:,2),A(:,1),A(:,3))
plot(A(:,2),A(:,3))