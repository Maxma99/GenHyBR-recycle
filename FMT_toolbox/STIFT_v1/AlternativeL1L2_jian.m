function [X] = AlternativeL1L2_jian(A,b,lambda1,N1,lambda2,N2)



L = speye(N*N);
options.lbound = 0;
options.ubound = 1;
k = 5;

for i = 1:k
    temp_x1 = l1F(A,b,lambda1,N1);
    b = b - A*temp_x1;
    temp_x2 = cglsTikConstraint(A,b,1:N2,lambda2,L,options);
    temp_x2 = temp_x2(:,end);
    b = b - A*temp_x2;
    X = X+temp_x1+temp_x2;
end
