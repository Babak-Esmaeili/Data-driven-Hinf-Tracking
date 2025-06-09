function [R,Q] = RQ_Gram(A)

if size(A,1) > size(A,2)
    A = A';
end
[m,n] = size(A);
% initialize Q and R
Q = zeros(n,m);
R = zeros(m,m);
% begin loop
for j = 1:m
    % calculate the values of R and Q for jth entry
    R(j,j) = sqrt(A(j,:)*A(j,:)');
    Q(:,j) = 1/R(j,j)*A(j,:)';
    % calculate the rest of the R values in the jth row
    R(j+1:m,j) = A(j+1:m,:)*Q(:,j);
    % update the A matrix
    A(j+1:m,:) = A(j+1:m,:) - R(j+1:m,j)*Q(:,j)';
end