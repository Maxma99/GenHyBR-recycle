function nrm = normM(v, M)
%
%     nrm = normM(v, M)
%
%  Calculate M norm for vector v.
%
% Input:
%          v - vector
%          M - matrix M for M norm
%
% Output:
%        nrm - result matrix
%
%  Reference:
%   Chung, de Sturler, and Jiang. "Hybrid Projection Methods with
%           Recycling for Inverse Problems". SISC, 2020.
%
% J. Chung, E. de Sturler, J. Jiang, and C. Ma, 2021



    if isa(M, 'function_handle')
      Mv = M(v);
    else
      Mv = M*v;
    end
    nrm = sqrt(v'*Mv);
    end