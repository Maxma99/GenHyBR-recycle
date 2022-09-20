function [x, IterInfo] = IRhybr(A, b, options)
%
%      [x, IterInfo] = IRhybr(A, b, options);
%
%  Hybrid Bidiagonalization Regularization for Iterative Image Restoration.
%  Note: This is a wrapper function for Julianne Chung's HyBR code.
%        See private/HyBR.m and associate codes, or 
%        http://www...
%
%   Input: A  -  object defining the coefficient matrix.
%          b  -  Right hand side vector.
% 
%   Optional Intputs:
%     options - Structure that can have:
%                x0      - initial guess; default is x0 = A'*b
%                x_true  - true solution, for testing purposes
%                MaxIter - integer specifying maximum number of iterations;
%                          default is min(size(A), 100)
%                Rtol    - stopping tolerance for the relative residual,
%                                norm(b - A*x)/norm(b)
%                          default is 1e-6
%
% 
%   Output:
%          x  -  solution
%   IterInfo  -  structure containing some information about the iteration
%                Iter     -  actual number of iterations performed
%                Rnrm     -  norm of the residual at each iteration
%                NE_Rnrm  -  scaled norm of the normal equations residual 
%                            at each iteration
%                Xnrm     -  norm of the solution at final iteration
%                Enrm     -  if input options specifies x_true, this contains
%                            the norm of the true relative error at final 
%                            iteration, norm(x - x_true)/norm(x_true)
%                StopFlag -  integer that indicates reason for stopping
%                            iteration:
%                               0 ==> Rtol satisfied
%                               1 ==> MaxIter reached
%
%                               
%
% J. Nagy, August, 2011
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%   [3] Chung, Nagy and O'Leary, "A Weighted-GCV Method for Lanczos-Hybrid
%       Regularization", ETNA 28 (2008), pp. 149-167.

% Initialization
defaultopt = struct('x0', 'none', 'MaxIter', [] ,'x_true', 'off', ...
  'Rtol', 1e-6, 'NE_Rtol', 1e-6, 'IterBar', 'on');

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 2
  error('IRhybr: Not Enough Inputs')
elseif nargin < 3
  options = [];
end
if isempty(options)
  options = defaultopt;
end

%
% Get options, and make sure everthing is initialized:
%
[m,n] = size(A);
defaultopt.MaxIter = min([m, n, 100]);
options = IRset(defaultopt, options);

x = IRget(options, 'x0', [], 'fast');
MaxIter = IRget(options, 'MaxIter', [], 'fast');
x_true = IRget(options, 'x_true', [], 'fast');
Rtol = IRget(options, 'Rtol', [], 'fast');
NE_Rtol = IRget(options, 'NE_Rtol', [], 'fast');
IterBar = IRget(options, 'IterBar', [], 'fast');

[mb, nb] = size(b);

%  
% For HyBR we'll make sure b to be a vector
if isa(A, 'single') || isa(A, 'double')
  if size(b(:),1) ~= m
    error('IRlhybr: Sizes of input matrix and vector are not compatible')
  else
      b = b(:);
  end
end
   
noIterBar = strcmp(IterBar,{'off'});

%
% SET HYBR OPTIONS and then extract from HyBRout the stuff 
% we need. 
HyBRoptions = HyBRset('x_true',x_true,'Iter',MaxIter);

[x, HyBRout] = modified_HyBR(A, b, [], HyBRoptions, noIterBar);

%x = reshape(x, mb, nb);

IterInfo.Iter = HyBRout.iterations;
IterInfo.Enrm = HyBRout.Enrm;
IterInfo.Rnrm = HyBRout.Rnrm;
IterInfo.Xnrm = HyBRout.Xnrm;
IterInfo.StopFlag = HyBRout.flag;





