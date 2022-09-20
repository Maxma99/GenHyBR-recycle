function [U, B, V, QV, H] = genrecyclingGKB(A, U, B, V, QV, Q, R, H, P_le, P_ri, W, Y, options)
%
%     [U, B, V, QV, H] = genrecyclingGKB(A, U, B, V, QV, Q, R, H, P_le, P_ri, W, Y, options)
%
%  Perform one step of generalized recycling GKB without reorthogonalization and 
%   without preconditioner.
%
% Input:
%          A - matrix
%       U, V - accumulation of vectors
%         QV - accumulation of Q*v_j
%          H - ADD definitions
%       Q, R - covariance matrices
%          B - bidiagonal matrix
% P_le, P_ri - inputs ignored
%       W, Y - Additional vectors for recycling
%    options - structure from HyBR (see HyBRset), ignored here
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%          H - updated upper triangular portion of the matrix
%
%  Reference:
%   Chung, de Sturler, and Jiang. "Hybrid Projection Methods with
%           Recycling for Inverse Problems". SISC, 2020.
%
% J. Chung, E. de Sturler, J. Jiang, and C. Ma, 2021

k = size(B,2)+1;

if k == 1
  Au = A'*(R^(-1)*U(:,k));
  WtAu = W'*Au;
  v = Au - W*WtAu;
else
  Au = A'*(R^(-1)*U(:,k));
  WtAu = W'*Au;  
  v = (Au - W*WtAu) - B(k, k-1)*V(:,k-1);
end
alpha = normM(v, Q);
v = v / alpha;
QV = [QV Q*v];
Av = A*QV(:,end); YtAv = Y'*(A*QV(:,end)); H(:,k) = YtAv;
u = Av - Y*YtAv - alpha*U(:,k);
beta = normM(u, @(x)R\x);
u = u / beta;
U = [U, u];
V = [V, v];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
