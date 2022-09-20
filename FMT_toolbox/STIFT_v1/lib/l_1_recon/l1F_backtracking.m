function [f] = l1F_backtracking(A,b,lambda,N)
%l1 regularization: find the solution to the problem 
%                  argmin_x{F(x):= ||Ax-b||^2+lambda*||x||_1}
%                  through the use of FISTA, a natural extension of the
%                  gradient-based method. Stepsize L_k for each iteraton is
%                  computed through backtracking.
%  Inputs:
%         eta: used to compute the stepsize
%         A: weighting matrix
%         b: measurements/response vector
%         lambda: regularization weight
%         N: number of iterations
% Output:
%         f: minimizer of the reularization problem
M=size(A,2);
x=zeros(M,1);
y=zeros(M,1);
t=1;
L=1;
eta=1.1;
ATA=A'*A;
ATb=A'*b;
for k=1:N
    %compute stepsie for step k
    i=0;
    l=eta^i*L;
    alpha=lambda/l;
    p=tau(alpha, y-2/l*(ATA*y-ATb));
    F=norm(A*p-b)^2+lambda*norm(p,1);
    Q=norm(A*y-b)^2+lambda*norm(p,1)+dot(p-y,2*(ATA*y-ATb));
    
    while F>Q
        i=i+1;
        l=eta^i*L;
        alpha=lambda/l;
        p=tau(alpha, y-2/l*(ATA*y-ATb));
        F=norm(A*p-b)^2+lambda*norm(p,1);
        Q=norm(A*y-b)^2+lambda*norm(p,1)+dot(p-y,2*(ATA*y-ATb));
    end
    %FISTA
    L=l;
    alpha=lambda/L;
    z=y-2/L*(ATA*y-ATb);
    xx=tau(alpha,z);
    tt=(1+sqrt(1+4*t^2))/2;
    y=xx+(t-1)/tt*(xx-x);
    x=xx;
    t=tt;
end
    f=y;
    
end

