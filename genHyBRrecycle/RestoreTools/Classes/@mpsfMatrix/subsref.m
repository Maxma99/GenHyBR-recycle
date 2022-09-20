function B = subsref(A, index)
%
%  Define field name indexing for mpsfMatrix object.
%
%               B = subsref(A, index)
%
%          This is called whenever a subscripted reference to the
%          psfMat2 object is made, such as:
%              A(i), i = 1, 2, 3, 4, 5, 6, 7
%              A.fieldname, fieldname = psf, matdata, type, boundary, transpose, center
%

%  J. Nagy 2/16/12

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = A.psf;
  case 2
    B = A.matdata;
  case 3
    B = A.type;
  case 4
    B = A.boundary;
  case 5
    B = A.transpose;
  case 6
    B = A.imsize;
  case 7
    B = A.center;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'psf'
    B = A.psf;
  case 'matdata'
    B = A.matdata;
  case 'type'
    B = A.type;
  case 'boundary'
    B = A.boundary;
  case 'transpose'
    B = A.transpose;
  case 'imsize'
    B = A.imsize;
  case 'center'
    B = A.center;
  otherwise
    B = subsref(A.psf, index);
  end

case '{}'
  error('Cell array indexing not supported for mpsfMatrix object.')
end
