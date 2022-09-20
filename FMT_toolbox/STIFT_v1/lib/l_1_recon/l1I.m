function [f] = l1I(A,b,lambda,N)
%l1 regularization: find the solution to the problem 
%                  argmin_x{F(x):= ||Ax-b||^2+lambda*||x||_1}
%                  through the use of ISTA, a natural extension of the
%                  gradient-based method
%  Inputs:
%         L: stepsize
%         A: design matrix
%         b: measurements/response vector
%         lambda: regularization weight
%         N: number of iterations
% Output:
%         f: minimizer of the reularization problem
M=size(A,2);
x=zeros(M,1);
L=2*eigs(A'*A,1);
alpha=lambda/L;
ATA=A'*A;
ATb=A'*b;
for k=1:N
    z=x-2/L*(ATA*x-ATb);
    x=tau(alpha,z);

end
    f=x;
end

