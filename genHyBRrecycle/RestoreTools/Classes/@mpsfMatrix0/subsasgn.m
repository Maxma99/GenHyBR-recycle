function A = subsasgn(A, index, val)
%
%  SUBSASGN  Define index assignment for mpsfMatrix0 object.
%
%               A = subsref(A, index, val)
%
%          This is called whenever an assignment statement of a
%          psf object is made, such as:
%              A(i) = val, i = 1, 2, 3, 4, 5, 6
%              A.fieldname = val,
%                fieldname = psf, matdata, transpose, imsize, 'center'
%

%  J. Nagy 2/23/12

switch index.type
case '()'
  switch index.subs{:}
  case 1
    A.psf = val;
  case 2
    A.matdata = val;
  case 5
    A.transpose = val;
  case 6
    A.imsize = val;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'psf'
    A.psf = val;
  case 'matdata'
    A.matdata = val;
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
  error('Cell array indexing not supported for mpsfMatrix0 object.')
end
