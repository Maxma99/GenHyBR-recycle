function A = mpsfMatrix0(varargin)
%
%  CONSTRUCTOR FOR mpsfMatrix0 (for multframe problems) OBJECT
%                  This only should be used for periodic boundary conditions.
%  
%  The mpsfMatrix0 class is based on a structure with four fields:
%    psf      - psf object
%    matdata  - matrix data needed to do matrix vector multiplications
%    transpose- indicates if the matrix has been transposed.
%    center   - location of point source
%
%  Calling Syntax:
%       A = mpsfMatrix0             (returns object with empty fields)
%       A = mpsfMatrix0(mpsfMatrix_Obj) (returns input object)
%       A = mpsfMatrix0( PSF )
%       A = mpsfMatrix0( PSF, center )
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
  A.transpose = 0;
  A.imsize = [];
  A.center = [];
  A = class(A, 'mpsfMatrix0');

case 1
  % if single argument of class mpsfMatrix0, return it
  % otherwise create an mpsfMatrix0
  % from a double array containing the PSF image
  %
  if ( isa( varargin{1}, 'psfMatrix0' ) )
    A = varargin{1};
  elseif ( isa(varargin{1}, 'double'))
    PSFs = varargin{1};
    A.psf = PSFs;
    center = round((size(PSFs(:,:,1))+1)/2);
    A.matdata = fft2(circshift(PSFs, 1-center));
    A.transpose = 0;
    A.center = center;
    A.imsize = [];

    A = class(A, 'mpsfMatrix0');
  else
    error('Incorrect argument type')
  end


case 2
  % if two argments, assume second input is the
  % center (a double array) 
  %
  PSFs = varargin{1};
  center = varargin{2};
  A.psf = PSFs;
  A.matdata = fft2(circshift(PSFs, 1-center));
  A.transpose = 0;
  A.center = center;
  A.imsize = [];

  A = class(A, 'mpsfMatrix0');

otherwise
  error('Incorrect number of input arguments.')
end


