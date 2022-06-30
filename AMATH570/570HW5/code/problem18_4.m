%This script verifies the theorem derived in part a of problem 18.4 in a
%specific case
clear
n = 5;
weights = ones(1,n+1); %given weights
A = zeros(n); 
A(1,2) = 1;
for k=1:n-2
    A(k+1,k)=k/(2*k+1); %note adjusted for indexing
    A(k+1,k+2)=(k+1)/(2*k+1);
end
A(n,n-1)=(n-1)/(2*n-1);

B = zeros(n);
B(n,:)=-n*weights(1:n)/((2*n-1)*weights(n+1)); %bottom row 

C = A+B;
fprintf('The eigenvalues of the comrade matrix are\n')
disp(eigs(C))

P = sum(legpoly(0:5),2);
fprintf('The roots of our polynomial are\n')
disp(roots(P,'all'))
