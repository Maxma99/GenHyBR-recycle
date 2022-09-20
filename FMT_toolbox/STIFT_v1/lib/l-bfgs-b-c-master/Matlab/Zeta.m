function [y] = Zeta(x)
% multi-dimensional version of scalar function zeta(x)=1/(x^2+delta_TV);
%   input a vector x of dimension n, output a vector of the same dimension
%   with zeta applied to each component of x.
% Remark: delta_TV is a predefined constant. 
delta_TV=10^-9;
Dim=size(x,1);
y=zeros(Dim,1);
for i=1:Dim
    y(i)=1/(x(i)^2+delta_TV);
end
end

