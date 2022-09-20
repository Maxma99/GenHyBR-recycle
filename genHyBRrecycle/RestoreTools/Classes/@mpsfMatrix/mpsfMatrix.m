function A = mpsfMatrix(varargin)
%
%  CONSTRUCTOR FOR mpsfMatrix (for multframe problems) OBJECT
%  
%  The mpsfMatrix class is based on a structure with four fields:
%    psf      - psf object
%    matdata  - matrix data needed to do matrix vector multiplications
%    type     - character string indicating:
%                 'invariant', 'variant', 'separable'
%    boundary - character string array indicating type of boundary conditions
%               to be used.  Choices are:
%                 'zero', 'periodic', 'reflexive' (or 'neumann'), 'antireflexive'
%               The default is reflexive.
%    transpose- indicates if the matrix has been transposed.
%
%  Calling Syntax:
%       A = mpsfMatrix             (returns object with empty fields)
%       A = mpsfMatrix(mpsfMatrix_Obj) (returns input object)
%       A = mpsfMatrix( PSF )
%       A = mpsfMatrix( PSF, center )
%       A = mpsfMatrix( PSF, boundary )
%       A = mpsfMatrix( PSF, center, boundary )
%       A = mpsfMatrix( PSF, boundary, center )
%
%    where 
%       * A   is an mpsfMatrix object
%       * PSF is a PSF
%       * b   is the right hand side image for the system Ax=b

%  J. Nagy  2/15/12

switch nargin

case 0
  A.psf = psf;
  A.matdata = [];
  A.type = '';
  A.boundary = '';
  A.transpose = 0;
  A.imsize = [];
  A.center = [];
  A = class(A, 'mpsfMatrix');

case 1
  % if single argument of class mpsfMatrix, return it
  % otherwise create an mpsfMatrix
  % from a double array containing the PSF image
  %
  if ( isa( varargin{1}, 'psfMatrix' ) )
    A = varargin{1};
  elseif ( isa(varargin{1}, 'double'))
    PSFs = varargin{1};
    A.psf = PSFs;
    Ppad = padarray(PSFs, [size(PSFs,1), size(PSFs,2), 0], 'post');
    center = round((size(PSFs(:,:,1))+1)/2);
    A.center = center;
    A.matdata = fft2(circshift(Ppad, 1-center));

    % Note that this should be used only for spatially invariant problems.
    A.type = 'invariant';
    A.boundary = 'reflexive';
    A.transpose = 0;
    A.imsize = [];

    A = class(A, 'mpsfMatrix');
  else
    error('Incorrect argument type')
  end


case 2
  % if two argments, need to check to see if the second one is the
  % center (a double array) or the boundary conditions (character string)
  %
  PSFs = varargin{1};
  if ( isa( varargin{2}, 'double') )
    boundary = 'reflexive';
    center = varargin{2};
  elseif ( isa( varargin{2}, 'char' ) )
    boundary = varargin{2};
    BC_types = {'zero', 'periodic', 'reflexive', 'neumann', 'antireflexive'};
    if ~any(strcmpi(BC_types, boundary))
      error(sprintf(['Incorrect boundary type: ', boundary]));
    end
    center = round((size(PSFs(:,:,1))+1)/2);
  end
  A.psf = PSFs;
  if strcmpi(boundary, 'periodic')
    Ppad = PSFs;
  else
    Ppad = padarray(PSFs, [size(PSFs,1), size(PSFs,2), 0], 'post');
  end
  A.matdata = fft2(circshift(Ppad, 1-center));

  % Note that this should be used only for spatially invariant problems.
  A.boundary = boundary;
  A.center = center;
  A.type = 'invariant';
  A.transpose = 0;
  A.imsize = [];

  A = class(A, 'mpsfMatrix');

case 3
  if ( isa( varargin{2}, 'double') )
    center = varargin{2};
    boundary = varargin{3};
  elseif ( isa( varargin{2}, 'char' ) )
    boundary = varargin{2};
    center = varargin{3};
  else
    error('incorrect input types')
  end
  BC_types = {'zero', 'periodic', 'reflexive', 'neumann', 'antireflexive'};
  if any(strcmpi(BC_types, boundary))
    A.boundary = boundary;
  else
    error(sprintf(['Incorrect boundary type: ', boundary]));
  end
  PSFs = varargin{1};
  A.psf = PSFs;
  if strcmpi(boundary, 'periodic')
    Ppad = PSFs;
  else
    Ppad = padarray(PSFs, [size(PSFs,1), size(PSFs,2), 0], 'post');
  end
  A.matdata = fft2(circshift(Ppad, 1-center));
  A.center =  center;
  % Note that this should be used only for spatially invariant problems.
  A.type = 'invariant';
  A.transpose = 0;
  A.imsize = [];

  A = class(A, 'mpsfMatrix');
       
otherwise
  error('Incorrect number of input arguments.')
end

