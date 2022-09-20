function [x,r] = L1LS(A,b,lambda,c)
%  2020-11-30  written by Jian
% Solve l1 regularized LS using Majorization-Minimization method, we set 
% c = 1e+5 in this problem. The Majorization-Minimization method for solving 
% problem updates x as:
%      x(k+1) =soft(1/c*A^T(y - Ax(k))+x(k),lambda/c),
%  inputs: 
%         A , b:  the equation
%         lambda: relax factor
%  outputs£º
%         x: the solution of LS

iter_time = 10000;
x = zeros(size(A,2),1);
% r = b;
for k = 1:iter_time
    s = x;
    temp = b-A*x;
%     r(k) =  norm(temp);
    x = soft(1/c.*A'*temp+x, lambda/c);
    if norm(x-s)<1e-5
        break;
    end   
end
end

%% soft is called the soft-thresholding operator
function z = soft(x, sigma)
    z = sign(x).*max(abs(x)- sigma,0);
end