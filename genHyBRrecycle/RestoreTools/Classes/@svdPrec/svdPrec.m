function P = svdPrec(varargin)
%
%  CONSTRUCTOR FOR svdPrec (svd preconditioner) OBJECT
%
%  The svdPrec class is based on a structure with four fields:
%    U - unitary matrix
%    S - vector containing singular or spectral values
%    V - unitary matrix
%
%  Calling Syntax:
%       P = svdPrec             (returns object with empty fields)
%       P = svdPrec(svdPrecObj) (returns input object)
%       P = svdPrec(  A, b )
%       P = svdPrec(  A, b, tol )
%       P = svdPrec(  A, b, 'help')
%
%    where 
%       * P   is a svdPrec object
%       * A   is a psfMatrix object or a multiPsfMatrix
%       * b   is the right hand side image for the system Ax=b that
%             is being preconditioned
%       * tol is a tolerance to "regularize" the preconditioner (e.g.,
%             the Hanke, Nagy, Plemmons approach)
%             If tol is not specified, a default will be chosen using
%             the generalized cross validation method.
%

%  J. Nagy  6/22/02

switch nargin

case 0
  P.u = [];
  P.s = [];
  P.v = [];
  P = class(P, 'svdPrec');

case 1
  if ( isa( varargin{1}, 'svdPrec' ) )
    P = varargin{1};
  else
    error('Incorrect argument type')
  end

case 2
  if ( isa( varargin{1}, 'psfMatrix') )
    [U, S, V] = svd(varargin{1}, varargin{2});
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'svdPrec');

  elseif (isa( varargin{1}, 'multiPsfMatrix') )
     error('svdPrec is not yet defined for multiPsfMatrix')

  else 
    error('Wrong input type')
  end    
        
case 3
  if ( isa( varargin{1}, 'psfMatrix') )
    [U, S, V] = svd(varargin{1}, varargin{2});
    if (ischar(varargin{3}))
      bhat = U'*varargin{2};
      trun_tol = GCVforSVD(S, bhat(:));
      trun_tol = trun_tol / max(S(:));
    else
      trun_tol = varargin{3};
    end
    S = S / max(S(:));
    S(S<trun_tol) = 1;
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'svdPrec');

  elseif ( isa ( varargin{1}, 'multiPsfMatrix') )
     error('svdPrec is not yet defined for multiPsfMatrix')
    
  else 
    error('Wrong input type')
  end   
  
otherwise
  error('Incorrect number of input arguments.')
end


