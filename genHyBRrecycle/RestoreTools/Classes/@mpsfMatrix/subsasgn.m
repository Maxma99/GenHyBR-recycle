function A = subsasgn(A, index, val)
%
%  SUBSASGN  Define index assignment for mpsfMatrix object.
%
%               A = subsref(A, index, val)
%
%          This is called whenever an assignment statement of a
%          psf object is made, such as:
%              A(i) = val, i = 1, 2, 3, 4, 5, 6, 7
%              A.fieldname = val,
%                fieldname = psf, matdata, type, boundary, transpose, imsize, center
%

%  J. Nagy 2/16/12

switch index.type
case '()'
  switch index.subs{:}
  case 1
    A.psf = val;
  case 2
    A.matdata = val;
  case 3
    A.type = val;
  case 4
    A.boundary = val;
  case 5
    A.transpose = val;
  case 6
    A.imsize = val;
  case 7
    A.center = val;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'psf'
    A.psf = val;
  case 'matdata'
    A.matdata = val;
  case 'type'
    A.type = val;
  case 'boundary'
    A.boundary = val;
  case 'transpose'
    A.transpose = val;
  case 'imsize'
    A.imsize = val;
  case 'center'
    A.center = val;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for mpsfMatrix object.')
end
