function [f] = l1_l2(A,b,lambda,N,gamma)
% joint l_1 and l_2 regularization: find the solution to the problem:
%                 argmin_x{F(x):= f(x)+g(x)}
%       f(x)=1/2*||Ax-b||^2; convex differentiable
%       g(x)=lambda*(||x||_1+gamma/2||beta||^2); convex nondifferentiable.
%       using Nesterov's accelerated method with constant stepsize.
% Inputs:
%         L: smallest Lipschitz constant of gradient(f), 
%            L(f)=lambda_max(A'A)
%         A: weighting matrix
%         b: measurements/response vector
%         lambda: regularization weight
%         N: number of iterations
%         gamma: weight of the ridge regression added to lasso
% Output:
%         f: minimizer of the reularization problem
M=size(A,2);
x=zeros(M,1);
y=zeros(M,1);
L=eigs(A'*A,1);
alpha=lambda/L;
beta=alpha*gamma;
ATA=A'*A;
ATb=A'*b;
for t=1:N
    z=y-1/L*(ATA*y-ATb);
    xx=1/(1+beta)*tau(alpha,z);
    y=xx+t/(t+3)*(xx-x);
    x=xx;
    
end
f=y;
end

